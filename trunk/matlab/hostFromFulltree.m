function host=hostFromFulltree(fulltree)
%Build vector 'host' indicating in which host each node of the fulltree is found
fathers=zeros(size(fulltree,1)+1,1);
fathers(fulltree(:,2)+1)=1:size(fulltree,1);
fathers(fulltree(:,3)+1)=1:size(fulltree,1);
fathers=fathers(2:end);
host=zeros(size(fulltree,1),1);
n=sum(fulltree(:,2)==0&fulltree(:,3)==0);
for i=1:n
    j=i;
    while 1
        host(j)=i;
        j=fathers(j);
        if fulltree(j,3)==0,break;end
    end
end
f=n+find(fulltree(n+1:end-1,3)==0);%All transmission events
for i=1:length(f)
    j=f(i);
    tocol=[];
    while host(j)==0
        tocol=[tocol j];
        j=fathers(j);
    end
    host(tocol)=host(j);
end