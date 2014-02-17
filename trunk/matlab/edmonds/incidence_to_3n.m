function[G]=incidence_to_3n(A);
% changes incidence matrix into 3 column format
N=size(A,1);
count=1;
for i=1:N
    for j=1:N
        if A(i,j)>0
        G(count,:)=[i,j,A(i,j)];
        count=count+1;
        end
        
    end
end




