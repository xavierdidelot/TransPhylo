%another example of usage
A=unidrnd(8,12,12)-1; % create a huge graph

for i=1:12
    A(i,i)=0; %no self loops allowed
end

V=[1:12]';
E=incidence_to_3n(A);

showgraph(E);
GT=edmonds(V,E);
reconstruct_2(GT)
