function [nodes,prob]=withinhost(times,neg)
%Input=times at which N samples are taken (counted forward in time from infection time)
%Output=array of size (2N)*3 where each row is a node, the first column indicate the date of the node and the last two
%columns indicate the two children
%This array has internal nodes sorted in order of most recent to most ancient node (and remains so during the algorithm)
%The last node corresponds to infection time and only has one child
%neg is the product of the within-host effective population size and the generation duration in days 
prob=0;
[tim,ind]=sort(times,'descend');
n=length(tim);
nodes=[0,ind(1),0];%Start with one node at time 0 and with the first isolate connected to it
i=2;
while i<=n%Graft branches one by one
    r=-log(rand)*neg;
    curt=tim(i);%Current time: start with date of isolate and go back in time until coalescence happens
    fi=find(nodes(:,1)<curt,1);
    for j=fi:size(nodes,1)
        if r>(curt-nodes(j,1))*(i-j)
            r=r-(curt-nodes(j,1))*(i-j);
            curt=nodes(j,1);
        else
            curt=curt-r/(i-j);%Found the time for grafting
            prob=prob+log(exppdf(r,neg));
            r=0;break;
        end
    end
    if r>0,continue;end%Rejection sampling forcing trees to be smaller than duration of infection
    %Create new node
    a=nodes(:,2:3);a(a>=j+n)=a(a>=j+n)+1;nodes(:,2:3)=a;%Renumbering according to table insertion in next line
    nodes=[nodes(1:j-1,:);curt ind(i) 0;nodes(j:end,:)];
    %Now choose on which branch to regraft amongst the branches alive at time curt
    no=j;side=2;
    prob=prob+log(1/(size(nodes,1)-j));
    w=1+floor(rand*(size(nodes,1)-j));
    while w>0
        no=no+side-1;side=3-side;
        if nodes(no,side+1)<=n || (nodes(no,side+1)>n && nodes(nodes(no,side+1)-n,1)>curt),w=w-1;end
    end
    nodes(j,3)=nodes(no,side+1);
    nodes(no,side+1)=n+j;
    i=i+1;
end
nodes=[zeros(n,3);nodes];
nodes(1:n,1)=times;