function fulltree=combine(ttree,tree)
%Combine a transmission tree and a phylogenetic tree into a fulltree

tree(:,1)=tree(:,1)-max(tree(:,1))+max(ttree(:,2));
n=ceil(size(tree,1)/2);
tree=[tree;zeros(n,3)];
source=find(ttree(:,3)==0);
tree(end,1)=ttree(source,1);
tree(end,2)=2*n-1;
notsource=setdiff(1:n,source);
i2=0;
for i=notsource
    i2=i2+1;
    f=find(tree(:,2)==i|tree(:,3)==i);
    fi=i;
    while tree(f,1)>ttree(i,1)
        fi=f;
        f=find(tree(:,2)==f|tree(:,3)==f);
        %if tree(f2,1)>ttree(i,1),fi=f;f=f2;else break;end 
    end
    if tree(f,2)==fi,tree(f,2)=2*n-1+i2;else tree(f,3)=2*n-1+i2;end
    tree(2*n-1+i2,2)=fi;
    tree(2*n-1+i2,1)=ttree(i,1);
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