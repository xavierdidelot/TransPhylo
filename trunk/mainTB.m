clear all
close all hidden force
seed=4;rand('state',seed);randn('state',seed);

%Parameters from which the MCMC is started
neg=100/365;%Within-host effective population size (Ne) times  generation duration (g)
popsize=100;%Size of population
beta=0.001;%Infectiveness
nu=1;%Rate of recovery

%Read TB data and tree
[id,month,day,year]=textread('tb/dates.csv','K%d\t%d/%d/%d');
dates=year+(month-1)/12+(day-1)/365;
timeLastRem=max(dates);
names=textread('tb/names.csv','%s');
ptree=phytreeread('tb/beast.nwk');
dist=get(ptree,'DISTANCES');
poin=get(ptree,'POINTERS');
nam=get(ptree,'LEAFNAMES');
for i=1:length(nam)
    names2{i}=names{str2num(nam{i})};
end
names=names2;
n=get(ptree,'NUMLEAVES');
ptree=zeros(n*2-1,3);
ptree(n+1:end,2:3)=poin;
for i=1:n*2-1
    time=0;
    j=i;
    while j<2*n-1
        time=time+dist(j);
        j=n+find(poin(:,1)==j|poin(:,2)==j);
    end
    ptree(i,1)=time;
end

if false
    load('matlabTB.mat');
else
%MCMC loop
fulltree=makeFullTreeFromPTree(ptree);%Starting point
mcmc=1000;
record(mcmc).tree=0;
pTTree=probTTree(ttreeFromFullTree(fulltree),popsize,beta,nu);
pPTree=probPTreeGivenTTree(fulltree,neg);
for i=1:mcmc
    if mod(i,100)==0,i,end
    %Record things
    record(i).tree=absolute(fulltree,timeLastRem);
    record(i).pTTree=pTTree;
    record(i).pPTree=pPTree;
    record(i).neg=neg;
    record(i).source=fulltree(end-1,4);
    record(i).nu=nu;
    record(i).beta=beta;
    %Metropolis update for transmission tree
    fulltree2=proposal(fulltree);
    pTTree2=probTTree(ttreeFromFullTree(fulltree2),popsize,beta,nu);
    pPTree2=probPTreeGivenTTree(fulltree2,neg);
    if log(rand)<pTTree2+pPTree2-pTTree-pPTree,fulltree=fulltree2;pTTree=pTTree2;pPTree=pPTree2;end
    %Metropolis update for Ne*g, assuming Exp(1) prior
    neg2=neg+(rand-0.5)*record(1).neg;
    if neg2<0,neg2=-neg2;end
    pPTree2=probPTreeGivenTTree(fulltree,neg2);
    if log(rand)<pPTree2-pPTree-neg2+neg,neg=neg2;pPTree=pPTree2;end
    %Gibbs updates for beta and nu based on O'Neill and Roberts JRSSA 16:121-129 (1999) and assuming Exp(1) priors
    ttree=ttreeFromFullTree(fulltree);
    intXtYtdt=intXtYtdtGivenTTree(ttree,popsize);
    beta=gamrnd(1+n-1,1/(1+intXtYtdt));
    nu=gamrnd(1+n,1/(1+n*mean(ttree(:,2)-ttree(:,1))));
end
end

ft=absolute(consensus(record,names),timeLastRem);
figure;
plotBothTree(ft,1,names);

%Plot traces
figure;
subplot(2,3,1);
plot([record(:).pPTree]+[record(:).pTTree]);
title('Log-posterior');
subplot(2,3,2);
plot([record(:).neg]);
title('Parameter Ne*g');
subplot(2,3,3);
plot([record(:).beta]);
title('Parameter beta');
subplot(2,3,4);
plot([record(:).nu]);
title('Parameter nu');
subplot(2,3,5);
plot([record(:).source]);
title('Source case');

%Posterior probability of being source case
subplot(2,3,6);
source=zeros(length(record)/2,1);
for i=1:length(source)
    source(i)=record(length(record)/2+i).tree(end-1,4);
end
hist([record(end/2:end).source],1:n);
title('Posterior probability of source case');