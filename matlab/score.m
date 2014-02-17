function s=score(record,truth)

n=sum(record(1).tree(:,2)==0&record(1).tree(:,3)==0);

%Analyze MCMC record
record=record(ceil(end/2):end);
mat=zeros(n+1,n+1);
for i=1:length(record)
    ttree=ttreeFromFullTree(record(i).tree);
    for j=1:size(ttree,1)
        mat(ttree(j,3)+1,j+1)=mat(ttree(j,3)+1,j+1)+1/length(record);
    end
end

%Analyze truth
matTruth=zeros(n+1,n+1);
ttree=ttreeFromFullTree(truth);
for j=1:size(ttree,1)
    matTruth(ttree(j,3)+1,j+1)=1;
end

%Calculate score
%s=mean(mat(find(matTruth==1)));
s=sum(mat(find(matTruth==1)))/sum(sum(mat));