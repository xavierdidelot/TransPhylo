% example of usage of the algorithms
V=[1:7]';
G=zeros(7);
G(1,2)=4;
G(1,3)=7;
G(1,5)=3;
G(1,7)=2;
G(2,6)=1;
G(3,4)=2;
G(3,5)=2;
G(4,6)=3;
G(5,4)=1;
G(6,5)=3;
G(4,2)=11;
G(7,6)=5;

 
E=incidence_to_3n(G);
GT= edmonds(V,E);

% Would need the bioinformatics toolbox for visualization
bg=biograph(G);bg.showWeights='On'; bg.view; 
% shows the maximum weight graph
TREEMAX=reconstruct_2(GT); 