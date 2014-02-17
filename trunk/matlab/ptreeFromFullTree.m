function ptree=ptreeFromFullTree(tree)
%Extract phylogenetic tree from a phylogenetic+transmission tree
n=sum(tree(:,2)+tree(:,3)==0);
tra=n+1;
while tra<size(tree,1)
    if tree(tra,3)~=0,tra=tra+1;continue;end
    t=tree(:,2:3);f=find(t==tra);t(f)=tree(tra,2);f=find(t>tra);t(f)=t(f)-1;tree(:,2:3)=t;
    tree=tree(setdiff(1:size(tree,1),tra),:);
    tra=n+1;
end
ptree=tree(1:end-1,1:end-1);