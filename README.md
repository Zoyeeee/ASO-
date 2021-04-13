# ASO-
题目:
---
![微信图片_20210413105538](https://user-images.githubusercontent.com/74950715/114490221-07a81880-9c47-11eb-8ac7-594a6b6352c4.jpg)
---
代码：
---
```matlab
clear all;clc
%%
%数据初始化
map_distance=[inf 3 1 5 8;3 inf 6 7 9;1 6 inf 4 2;5 7 4 inf 3;8 9 2 3 inf];
Ng=100;%迭代次数
antnum=30;%蚂蚁数目
Alpha=1;%信息素重要程度[0,5]
Beta=2;%启发因子重要程度[0,5]
Rho=0.75;%信息素蒸发[0.1,0.99]
Q=1;%信息素增加强度
Eta=1./map_distance;%启发因子
%%
%初始化蚁群
ctitynum=size(map_distance,1);%城市数量
BestRoute=nan(1,ctitynum);
BestDistence=inf;
RouteOfAnt=nan(antnum,ctitynum);%记录蚂蚁的路线
for i=1:antnum
    RouteOfAnt(i,:)=randperm(ctitynum);%初始化每只蚂蚁的位置
end
PheromoneMap=zeros(ctitynum,ctitynum);%初始化信息素矩阵
%%
for i=1:Ng
    [RouteOfAnt,PheromoneMap,BestRoute,BestDistence]=SelectProbability(Alpha,Beta,RouteOfAnt,PheromoneMap,Eta,map_distance,BestDistence,BestRoute,Rho,Q);
end
BestRoute
BestDistence
%%
function [RouteOfAnt,PheromoneMap,BestRoute,BestDistence]=SelectProbability(Alpha,Beta,RouteOfAnt,PheromoneMap,Eta,map_distance,BestDistence,BestRoute,Rho,Q)
    ctitynum=size(PheromoneMap,1);
    antnum=size(RouteOfAnt,1);
    Probability=nan(1,ctitynum);
    trueSelect=[];
    for j=1:ctitynum
        for i=1:antnum
            if j~=ctitynum
                for k=1:ctitynum
                    if ismember(k,RouteOfAnt(i,:))
                        Probability(k)=0;
                    else
                        Probability(k)=(Eta(j,k)^Beta)*(PheromoneMap(j,k)^Alpha);
                    end
                end
                if sum(Probability)~=0
                    Probability=Probability/(sum(Probability));
                    Pcum=cumsum(Probability);
                    Select=find(Pcum>=rand());
                    if size(Select)~=1
                        trueSelect=Select(end);
                    else
                        trueSelect=Select;
                    end
                    RouteOfAnt(i,j+1)=trueSelect;
                end
            else
                [BestRoute,BestDistence]=DistanceComputer(RouteOfAnt,map_distance,BestDistence,BestRoute);%
                [PheromoneMap]=PheromoneMapUpdate(PheromoneMap,Rho,antnum,RouteOfAnt,Q,map_distance);
            end 
        end
    end
end
function [BestRoute,BestDistence]=DistanceComputer(RouteOfAnt,map_distance,BestDistence,BestRoute)
    antnum=size(RouteOfAnt,1);
    ctitynum=size(RouteOfAnt,2);
    Distance_save=nan(1,antnum);
    tempsum=0;
    for i=1:antnum
        for j= 1:(ctitynum-1)
            tempsum=tempsum+map_distance(RouteOfAnt(i,j),RouteOfAnt(i,j+1));
        end
        Distance_save(1,i)=tempsum;
    end
    tempmin=min(Distance_save);
    tempminID=find(Distance_save==tempmin);
    if tempmin<BestDistence
        BestDistence=tempmin;
        BestRoute=RouteOfAnt(tempminID,:);
    end
end
function [PheromoneMap]=PheromoneMapUpdate(PheromoneMap,Rho,antnum,RouteOfAnt,Q,map_distance)
    ctitynum=size(PheromoneMap,1);
    Delta_Tau = zeros(ctitynum, ctitynum);
    for i = 1: antnum
        for j = 1: (ctitynum - 1)
            Delta_Tau(RouteOfAnt(i, j), RouteOfAnt(i, j + 1)) = Delta_Tau(RouteOfAnt(i, j), RouteOfAnt(i, j + 1)) + Q / map_distance(RouteOfAnt(i, j), RouteOfAnt(i, j + 1));
        end
        Delta_Tau(RouteOfAnt(i, 1), RouteOfAnt(i, ctitynum)) = Delta_Tau(RouteOfAnt(i, 1), RouteOfAnt(i, ctitynum)) + Q / map_distance(RouteOfAnt(i, j), RouteOfAnt(i, j + 1));
    end
    PheromoneMap = (1 - Rho) .* PheromoneMap + Delta_Tau;
end
```
结果:
---
运行多次后得到：
路线：2-1-3-5-4
9
