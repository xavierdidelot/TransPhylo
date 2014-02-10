function[]=showgraph(G)
% takes input as structure from edmonds
if size(G,2)==3 
    N=max(unique(G(:,1:2)));
    G=[G;N N 0];
INCIDENCE=sparse(G(:,1),G(:,2),G(:,3));
elseif size(G,1)==size(G,2)
    INCIDENCE=sparse(G); 
end

bg=biograph(INCIDENCE);bg.showWeights='On'; 
bg.view;
