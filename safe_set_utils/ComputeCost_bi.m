function [ Qfun] = ComputeCost_bi( x,u,Q)
s_f=40;

for i = 1:(size(x,2)-1)
    Index = size(x,2)-1-i+1; % Need to start from the end
    if i == 1
        Cost(Index) = Q*max(s_f-x(1,Index),0);    
    else
        Cost(Index) = Cost(Index+1) + Q*max(s_f-x(1,Index),0);    
    end
end
Qfun = Cost;
end