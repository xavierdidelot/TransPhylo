function s=scoreConsensus(cons,truth)

n=sum(truth(:,2)==0&truth(:,3)==0);

%Analyze consensus
mat=zeros(n+1,n+1);
ttree=ttreeFromFullTree(cons);
for j=1:size(ttree,1)
    mat(ttree(j,3)+1,j+1)=1;
end

%Analyze truth
matTruth=zeros(n+1,n+1);
ttree=ttreeFromFullTree(truth);
for j=1:size(ttree,1)
    matTruth(ttree(j,3)+1,j+1)=1;
end

%Calculate score
s=mean(mat(find(matTruth==1)));