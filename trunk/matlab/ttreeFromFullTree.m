function ttree=ttreeFromFullTree(fulltree)
%Takes in a fulltree, and extracts the ttree
host=fulltree(:,4);
ttree=fulltree(fulltree(:,2)==0&fulltree(:,3)==0,1);
n=length(ttree);
ttree=[zeros(length(ttree),1) ttree zeros(length(ttree),1)];
for i=1:n
    j=i;
    while host(j)==i
        j=find(fulltree(:,2)==j | fulltree(:,3)==j);
    end
    ttree(i,1)=fulltree(j,1);
    ttree(i,3)=host(j);
end
