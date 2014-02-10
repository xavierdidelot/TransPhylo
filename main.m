clear all
close all hidden force
rand('state',0);
randn('state',0);

%Biological parameters, calibrated for Staphylococcus aureus,  based on 
%Young et al. (2012) PNAS 109:4550-4555
%Golubchik et al. (2013) PLoS One 8:e61319
neg=100/365;%Within-host effective population size (Ne) times  generation duration (g)
popsize=1000;%Size of population
beta=0.001;%Infectiveness
nu=1;%Rate of recovery

%Create a transmission tree with at least ten individuals
n=1;
while n~=5%n<10 || n>15
    n
    ttree=makeTTree(popsize,beta,nu);
    n=size(ttree,1);
end

%Create a within-host phylogenetic tree for each infected host
for i=1:n
    times=[ttree(i,2);ttree(find(ttree(:,3)==i),1)]-ttree(i,1);
    wtree{i}=withinhost(times,neg);
end

%Glue these trees together
truth=glueTrees(ttree,wtree);
truth(:,1)=truth(:,1)+2005;%Epidemic started in 2005
timeLastRem=max(truth(:,1));

return;

%plotBothTree(truth,1);
%return;

%MCMC loop
fulltree=makeFullTreeFromPTree(ptreeFromFullTree(truth));%Starting point
mcmc=10000;
record(mcmc).tree=0;
pTTree=probTTree(ttreeFromFullTree(fulltree),popsize,beta,nu);
pPTree=probPTreeGivenTTree(fulltree,neg);
for i=1:mcmc
    i
    %if sanitycheck(fulltree)==0,break;end
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
    %%Metropolis update for Ne*g, assuming Exp(1) prior
    %neg2=neg+(rand-0.5)*record(1).neg;
    %if neg2<0,neg2=-neg2;end
    %pPTree2=probPTreeGivenTTree(fulltree,neg2);
    %if log(rand)<pPTree2-pPTree-neg2+neg,neg=neg2;pPTree=pPTree2;end
    %%Gibbs updates for beta and nu based on O'Neill and Roberts JRSSA 16:121-129 (1999) and assuming Exp(1) priors
    %ttree=ttreeFromFullTree(fulltree);
    %intXtYtdt=intXtYtdtGivenTTree(ttree,popsize);
    %beta=gamrnd(1+n-1,1/(1+intXtYtdt));
    %nu=gamrnd(1+n,1/(1+n*mean(ttree(:,2)-ttree(:,1))));
end

ft=absolute(consensus(record),timeLastRem);

%Plot four trees
figure;
subplot('Position',[0.1 0.5 0.35 0.35]);
plotTTree(ttreeFromFullTree(truth));
%title('Correct transmission tree','FontSize',16);
subplot('Position',[0.5 0.5 0.35 0.35]);
set(gca,'FontSize',12);
plotBothTree(truth,1);
%title('Correct coloring','FontSize',16);
subplot('Position',[0.1 0.1 0.35 0.35]);
plotTTree(ttreeFromFullTree(ft));
%title('Inferred transmission tree','FontSize',16);
subplot('Position',[0.5 0.1 0.35 0.35]);
set(gca,'FontSize',12);
plotBothTree(ft,1);
%title('Inferred coloring','FontSize',16);
%return;

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