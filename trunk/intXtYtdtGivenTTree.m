function intXtYtdt=intXtYtdtGivenTTree(ttree,popsize)
%return the integral of Xt*Yt*dt

[times,ind]=sort([ttree(:,1);ttree(:,2)]);%Event times
types=ind>size(ttree,1);%0 for infection, 1 for recovery
sir=[popsize-1 1 0];
intXtYtdt=0;
for i=2:length(times)
    intXtYtdt=intXtYtdt+sir(1)*sir(2)*(times(i)-times(i-1));
    if types(i)==1
        %Recovery
        sir(3)=sir(3)+1;sir(2)=sir(2)-1;
    else
        %Infection
        sir(1)=sir(1)-1;sir(2)=sir(2)+1;
    end
end