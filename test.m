clear all
close all hidden force
rand('state',0);
randn('state',0);

neg=1/365;%Within-host effective population size (Ne) times  generation duration (g)
R=1;%Basic reproduction number

%Create a transmission tree with ten individuals
n=1;
while n~=10
    [ttree,p1]=makeTTree(R);
    n=size(ttree,1);
end

%Create a within-host phylogenetic tree for each infected host
for i=1:n
    times=[ttree(i,2);ttree(find(ttree(:,3)==i),1)]-ttree(i,1);
    [wtree{i},p]=withinhost(times,neg);
    p1=p1+p;
end


%Glue these trees together
truth=glueTrees(ttree,wtree);
%truth(:,1)=truth(:,1)+2005;%Epidemic started in 2005
timeLastRem=max(truth(:,1));
[p1 probTTree(ttreeFromFullTree(truth),R)+probPTreeGivenTTree(truth,neg)]
