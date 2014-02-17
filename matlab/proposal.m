function res=proposal(fulltree)
n=sum(fulltree(:,2)==0&fulltree(:,3)==0);
host=fulltree(:,4);
fathers=zeros(size(fulltree,1)+1,1);
fathers(fulltree(:,2)+1)=1:size(fulltree,1);
fathers(fulltree(:,3)+1)=1:size(fulltree,1);
fathers=fathers(2:end);

%Choose a transmission event
f=find(fulltree(:,2)>0&fulltree(:,3)==0);
w=1+floor(rand*length(f));
w=f(w);

if w==size(fulltree,1)
    %Move for the transmission to the index case
    r=(rand-0.5)/1;
    mini=-fulltree(end-1,1);
    if r<mini,r=mini+(mini-r);end
    fulltree(1:end-1,1)=fulltree(1:end-1,1)+r;
else
    %Move for the transmission to someone else
    infector=host(w);
    infected=host(fulltree(w,2));
    %First remove the transmission node
    f=fathers(w);
    if fulltree(f,2)==w,fulltree(f,2)=fulltree(w,2);fathers(fulltree(w,2))=f;else fulltree(f,3)=fulltree(w,2);fathers(fulltree(w,2))=f;end
    fulltree(w,:)=0;
    %Second add it back on the path from infector to infected
    islocs=zeros(size(fulltree,1),1);
    path1=infector;islocs(path1)=1;
    while path1<size(fulltree,1),path1=fathers(path1);islocs(path1)=1;end;
    path2=infected;islocs(path2)=1;
    while path2<size(fulltree,1),path2=fathers(path2);islocs(path2)=1-islocs(path2);end
    locs=find(islocs);
    lens=zeros(length(locs),1);
    for i=1:length(locs)
        lens(i)=fulltree(locs(i),1)-fulltree(fathers(locs(i)),1);
    end
    r=rand*sum(lens);
    for i=1:length(locs)
        if r>lens(i)
            r=r-lens(i);
        else
            fulltree(w,1)=fulltree(locs(i),1)-r;
            j=fathers(locs(i));
            if fulltree(j,2)==locs(i),fulltree(j,2)=w;else fulltree(j,3)=w;end
            fulltree(w,2)=locs(i);            
            break;
        end
    end
end

%Reorder nodes chronologically
[~,ind]=sort(fulltree(n+1:end,1),'descend');
for i=n+1:size(fulltree,1)
    for j=2:3
        if fulltree(i,j)>n,fulltree(i,j)=n+find(ind==fulltree(i,j)-n);end
    end
end
fulltree=fulltree([(1:n)';n+ind],:);
res=[fulltree(:,1:3) hostFromFulltree(fulltree)];