function prob=probTTree(ttree,popsize,beta,nu)
%Calculate the log-probability of a transmission tree given popsize, beta and nu

for i=1:size(ttree,1)
    infector=ttree(i,3);
    if infector==0,continue;end
    if ttree(i,1)>ttree(infector,2)||ttree(i,1)<ttree(infector,1),prob=-inf;return;end
end

[times,ind]=sort([ttree(:,1);ttree(:,2)]);%Event times
types=ind>size(ttree,1);%0 for infection, 1 for recovery
sir=[popsize-1 1 0];
prob=0;
for i=2:length(times)
    rate=nu*sir(2)+beta*sir(1)*sir(2);
    prob=prob+log(rate)-rate*(times(i)-times(i-1));%;log(exppdf(times(i)-times(i-1),1/rate));
    prob=prob+log(1/sir(2));
    if types(i)==1
        %Recovery
        prob=prob+log(nu*sir(2)/rate);
        sir(3)=sir(3)+1;sir(2)=sir(2)-1;
    else
        %Infection
        prob=prob+log(1-nu*sir(2)/rate);
        sir(1)=sir(1)-1;sir(2)=sir(2)+1;
    end
end