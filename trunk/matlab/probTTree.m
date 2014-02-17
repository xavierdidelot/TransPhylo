function prob=probTTree(ttree,R)
%Calculate the log-probability of a transmission tree given R
%Assume w is Gamma(2,1)
prob=0;
for i=1:size(ttree,1)
    prob=prob+log(gampdf(ttree(i,2)-ttree(i,1),2,1));
    offspring=find(ttree(:,3)==i);
    prob=prob+log(poisspdf(length(offspring),R));
    for j=offspring'
        prob=prob+log(gampdf(ttree(j,1)-ttree(i,1),2,1));
    end
end