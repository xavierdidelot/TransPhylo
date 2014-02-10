function [ttree,prob]=makeTTree(popsize,beta,nu)
%Creates a transmission tree, and returns the result as a N*3 matrix in the
%following format:
%One row per infected host
%First column is time of infection
%Second column is time of recovery
%Third column is infector
%nu is the rate of recovery
%beta is the rate of infection
%popsize is the total size of population
ttree=[0 0 0];
sir=[popsize-1 1 0];
curt=0;
prob=0;
while 1
    rate=nu*sir(2)+beta*sir(1)*sir(2);
    if rate==0,break;end
    t=exprnd(1/rate);
    prob=prob+log(exppdf(t,1/rate));
    curt=curt+t;
    w=find(ttree(:,2)==0);%All infected
    prob=prob+log(1/length(w));
    w=w(1+floor(length(w)*rand));%Pick one at random
    if rand<nu*sir(2)/rate
        %Recovery
        prob=prob+log(nu*sir(2)/rate);
        ttree(w,2)=curt;
        sir(3)=sir(3)+1;sir(2)=sir(2)-1;
    else
        %Infection
        prob=prob+log(1-nu*sir(2)/rate);
        ttree=[ttree;curt 0 w];
        sir(1)=sir(1)-1;sir(2)=sir(2)+1;
    end
end