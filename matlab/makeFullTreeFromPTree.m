function fulltree=makeFullTreeFromPTree(tree)
%Return a non-zero probability phylogenetic+transmission tree for a given phylogenetic tree
n=ceil(size(tree,1)/2);
tree=[tree;zeros(n,3)];
tree(end,1)=min(tree(1:(2*n-1),1))-1;
tree(end,2)=2*n-1;
[~,source]=max(tree(1:n,1));%The source is the last leaf
notsource=setdiff(1:n,source);
i2=0;
for i=notsource
    i2=i2+1;
    f=find(tree(:,2)==i|tree(:,3)==i);
    if tree(f,2)==i,tree(f,2)=2*n-1+i2;else tree(f,3)=2*n-1+i2;end
    tree(2*n-1+i2,2)=i;
    tree(2*n-1+i2,1)=(tree(f,1)+tree(i,1))/2;
end
tree(:,1)=tree(:,1)-min(tree(:,1));

%Reorder nodes chronologically
[~,ind]=sort(tree(n+1:end,1),'descend');
for i=n+1:size(tree,1)
    for j=2:3
        if tree(i,j)>n,tree(i,j)=n+find(ind==tree(i,j)-n);end
    end
end
tree=tree([(1:n)';n+ind],:);
fulltree=[tree hostFromFulltree(tree)];