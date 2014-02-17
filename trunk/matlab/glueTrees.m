function fulltree=glueTrees(ttree,wtree)
%Glue together the within-host trees using the transmission tree
n=size(ttree,1);
no=n+1;
for i=n:-1:1%Consider all hosts
    ni=size(wtree{i},1)/2;
    for j=1:ni
        lab(i,j)=no;no=no+1;
    end
    labh(i)=no-1;
end
leaves=[];
intnodes=[];
for i=1:n
    tree=wtree{i};
    ni=size(wtree{i},1)/2;
    tree(:,1)=tree(:,1)+ttree(i,1);%Add infection time to all nodes
    leaves=[leaves;tree(1,:)];
    tree=tree(ni+1:end,:);
    f=find(ttree(:,3)==i);%Infected by current nodes
    for j=1:size(tree,1)
        for k=2:3
            if tree(j,k)==0,continue;end
            if tree(j,k)==1,tree(j,k)=i;
            elseif tree(j,k)<=ni,tree(j,k)=labh(f(tree(j,k)-1));
            else tree(j,k)=lab(i,tree(j,k)-ni);end
        end
    end
    intnodes=[tree;intnodes];
end
fulltree=[leaves;intnodes];
fulltree=[fulltree hostFromFulltree(fulltree)];
