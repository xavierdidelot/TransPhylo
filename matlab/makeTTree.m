function [ttree,prob]=makeTTree(R)
%Creates a transmission tree, and returns the result as a N*3 matrix in the
%following format:
%One row per infected host
%First column is time of infection
%Second column is time of sampling
%Third column is infector
%R is the basic reproduction number
%Assume w is Gamma(2,1)
ttree=[0 0 0];
prob=0;
todo=1;
while ~isempty(todo)
    draw=gamrnd(2,1);
    ttree(todo(1),2)=ttree(todo(1),1)+draw;
    prob=prob+log(gampdf(draw,2,1));
    offspring=poissrnd(R);
    prob=prob+log(poisspdf(offspring,R));
    for i=1:offspring
        draw=gamrnd(2,1);
        prob=prob+log(gampdf(draw,2,1));
        ttree=[ttree;ttree(todo(1),1)+draw 0 todo(1)];
        if size(ttree,1)>100,ttree=[];return;end
        todo=[todo size(ttree,1)];
    end
    todo=todo(2:end);
end