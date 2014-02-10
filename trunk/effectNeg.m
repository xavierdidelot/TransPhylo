clear all
close all hidden force
seed=89;
rand('state',seed);
randn('state',seed);

neg=100/365;%Within-host effective population size (Ne) times  generation duration (g)
popsize=100;%Size of population
beta=0.02;%Infectiveness
nu=2;%Rate of recovery

%Create a transmission tree with at least ten individuals
n=1;
while n~=10
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

neg=0.001;
%Create a within-host phylogenetic tree for each infected host
for i=1:n
    times=[ttree(i,2);ttree(find(ttree(:,3)==i),1)]-ttree(i,1);
    wtree{i}=withinhost(times,neg);
end

%Glue these trees together
truth2=glueTrees(ttree,wtree);
truth2(:,1)=truth2(:,1)+2005;%Epidemic started in 2005

%Plot three trees
figure;
subplot('Position',[0.1 0.6 0.8 0.4]);
plotTTree(ttreeFromFullTree(truth));
subplot('Position',[0.1 0.4 0.8 0.2]);
xlim([2005 2008]);set(gca,'FontSize',12);set(gca,'XTick',[2005 2006 2007 2008 ]);
%xs=xlim;ys=ylim;text(xs(1),ys(2),'A');
plotBothTree(truth2,1);
subplot('Position',[0.1 0.1 0.8 0.2]);
xlim([2005 2008]);set(gca,'FontSize',12);set(gca,'XTick',[2005 2006 2007 2008 ]);
plotBothTree(truth,1);
   ax = axes('position',[0,0,1,1],'visible','off');
   tx = text(0.1,0.95,'A','FontSize',20);
   tx = text(0.1,0.6 ,'B','FontSize',20);
   tx = text(0.1,0.3 ,'C','FontSize',20);

