function plotPTree(tree,method)
%This function takes a fulltree as input
if nargin==1,method=2;end
n=sum(tree(:,2)+tree(:,3)==0);
a=gca;

if method==1
    %Add ghost leaves to represent transmission nodes
    tra=[];
    for i=n+1:size(tree,1)-1
        if tree(i,3)==0,tra=[tra i];end
    end
    t=tree(:,2:3);t(t>n)=t(t>n)+length(tra);tree(:,2:3)=t;
    no=n+1;
    toadd=[];
    for i=n+1:size(tree,1)-1
        if tree(i,3)==0,tree(i,3)=no;no=no+1;toadd=[toadd;tree(i,1) 0 0 0];end
    end
    tree=[tree(1:n,:);toadd;tree(n+1:end,:)];
else
    %Remove transmission nodes
    tree=ptreeFromFullTree(tree);
end

n=ceil(size(tree,1)/2);
br=tree(n+1:end,2:3);
for i=1:n*2-2
    f=find(tree(:,2)==i|tree(:,3)==i);
    d(i)=tree(i,1)-tree(f,1);
end
d=[d 0];
plot(phytree(br,d'));
copyobj(allchild(gca),a);
close;
axes(a);
axis off;