function prob=probPTreeGivenTTree(fulltree,neg)
%To calculate the probability of a phylogenetic tree given a transmission
%tree, split up the phylogenetic tree into individual within-host subtrees, and
%calculate the product of their probabilities
n=sum(fulltree(:,2)==0&fulltree(:,3)==0);
prob=0;
fathers=zeros(size(fulltree,1)+1,1);
fathers(fulltree(:,2)+1)=1:size(fulltree,1);
fathers(fulltree(:,3)+1)=1:size(fulltree,1);
fathers=fathers(2:end);
for i=1:n
    subtree=extractSubtree(fulltree,i,fathers);
    prob=prob+probSubtree(subtree,neg);
end

function subtree=extractSubtree(fulltree,which,fathers)
%Take all nodes in host
host=fulltree(:,4);
ind=1:size(fulltree,1);
ind=ind(host==which);
%Add father of oldest node
[~,j]=min(fulltree(ind));
ind=[ind fathers(ind(j))];
%Create subtree
subtree=zeros(length(ind),2);
invind=zeros(size(fulltree,1));invind(ind)=1:length(ind);
for i=1:length(ind)
    subtree(i,1)=fulltree(ind(i),1);
    if i<length(ind)
        subtree(i,2)=invind(fathers(ind(i)));
    end
end

function p=probSubtree(tab,rate)
%tab(:,1)=age at bottom;tab(:,2)=father;rate=coalescence rate
%Return the log-prior probability of a subtree
%This is an extension to Eq1 of Drummond et al (2002) Genetics
%161:1307-1320 that accounts for condition TMRCA<INCUBATION_PERIOD
p=0;
tab(:,1)=max(tab(:,1))-tab(:,1);%convert times to ages
isiso=ones(size(tab,1),1);
isiso(tab(1:end-1,2))=0;
iso=find(isiso);
ex=zeros(size(tab,1),1);%Keep track of which nodes are active
[s,ind]=sort(tab(:,1));
cur=iso(1);while tab(cur,2)>0,ex(cur)=1;cur=tab(cur,2);end;ex(end)=1;%Activate all ancestors of first leave
for l=2:length(iso)%For all others leaves in increasing order of age
    isanc=zeros(size(tab,1),1);%Ancestors of the current node
    anc=iso(l);while tab(anc,2)>0,isanc(anc)=1;anc=tab(anc,2);end
    bra1=0;
    bra2=0;
    start=false;
    found=false;
    curage=0;
    k=0;
    for i=1:length(s)
        if ind(i)==iso(l),start=true;curage=tab(iso(l),1);end
        if ex(ind(i))==0,continue;end%Ignore non-existant nodes
        if start
            if found==false
                bra1=bra1+k*(tab(ind(i),1)-curage);
                if isanc(ind(i)),found=true;end
            end
            bra2=bra2+k*(tab(ind(i),1)-curage);
        end
        curage=tab(ind(i),1);
        if isiso(ind(i)),k=k+1;else k=k-ex(ind(i))+1;end
    end
    p=p-log(rate)-bra1/rate-log(1-exp(-bra2/rate));
    if k~=1,'errorHere',end
    if bra1>=bra2,'errorThere',end
    cur=iso(l);while tab(cur,2)>0,if ex(cur)==1,ex(cur)=2;break;end;ex(cur)=1;cur=tab(cur,2);end%Make all ancestors of current node active
end

function p=probSubtreeOLD(tab,rate)
%Return the log-prior in Eq1 of Drummond et al (2002) Genetics 161:1307-1320
%This is not quite right because of the truncation of the tree distribution TMRCA<INCUBATION_PERIOD
n=size(tab,1)/2;
p=-log(rate)*(n-1);
tab(:,1)=max(tab(:,1))-tab(:,1);%convert times to ages
[s,ind]=sort(tab(:,1));
k=1;
for i=2:length(s)-1
    p=p-(k*(k-1)/(2*rate)*(s(i)-s(i-1)));
    if isempty(find(tab(:,2)==ind(i))),k=k+1;else k=k-sum(tab(:,2)==ind(i))+1;end
end
if k~=1,'errorP',end