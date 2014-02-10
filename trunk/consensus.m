function fulltree=consensus(record,names)

n=sum(record(1).tree(:,2)==0&record(1).tree(:,3)==0);
if nargin<2,for i=1:n,names{i}=sprintf('%d',i);end;end
plotting=false;

%Analyze MCMC record
record=record(ceil(end/2):end);
mat=zeros(n+1,n+1);
tinf=zeros(n,1);
trec=zeros(n,1);
tinf2=zeros(n+1,n+1,2);
for i=1:length(record)
    ttree=ttreeFromFullTree(record(i).tree);
    for j=1:size(ttree,1)
        tinf(j)=tinf(j)+ttree(j,1)/length(record);
        trec(j)=trec(j)+ttree(j,2)/length(record);
        mat(ttree(j,3)+1,j+1)=mat(ttree(j,3)+1,j+1)+1/length(record);
        tinf2(ttree(j,3)+1,j+1,1)=tinf2(ttree(j,3)+1,j+1,1)+ttree(j,1);
        tinf2(ttree(j,3)+1,j+1,2)=tinf2(ttree(j,3)+1,j+1,2)+1;
    end
end

%Plot graph
mat2=mat;
if plotting==true
ids{1}='SOURCE';
for i=2:n+1
    ids{i}=sprintf('%s [%.2f;%.2f]',names{i-1},tinf(i-1)-000,trec(i-1)-000);
end
mat(mat<0.1)=0;
bg=biograph(mat,ids,'LayoutScale',0.3);
colorEdges(bg,ids,mat2);
g=biograph.bggui(bg);
fig=figure;a=gca;
copyobj(allchild(g.biograph.hgAxes),a);%Copy biograph viewer to current axes
set(0,'ShowHiddenHandles','on');close(setdiff(get(0,'Children'),fig));%Close the biograph viewer
axes(a);
axis off
set(gcf,'Color','w');
end

%Find best tree using Edmonds' algorithm from Edmonds (1967) J. Res. Nat. Bur. Standards 71B:233-240
cd edmonds
[~,ind]=max(mat2(1,:));val=mat2(1,ind);mat2(1,setdiff(1:end,ind))=0;mat2(1,ind)=1;%Only keep most likely source, otherwise may be more than one child of the source
E=incidence_to_3n(mat2);
mat2(1,ind)=val;
GT= edmonds([1:n+1]',E);
[~,GTM]=reconstruct_2(GT);
cd ..
ids{1}='SOURCE';
ttree=zeros(n,3);
for i=2:n+1
    fi=find(GTM(:,i)>0);%Infector
    if length(fi)>1,sprintf('More than one infector for %d!',i);end
    if isempty(fi),sprintf('No infector for %d! %f',i,sum(mat2(:,i+1))),continue;end
    ids{i}=sprintf('%s [%.1f;%.1f]',names{i-1},tinf2(fi,i,1)/tinf2(fi,i,2)-000,trec(i-1)-000);
    ttree(i-1,1)=tinf2(fi,i,1)/tinf2(fi,i,2);
    ttree(i-1,2)=trec(i-1);
    ttree(i-1,3)=fi-1;
end

%Plot tree
if plotting==true
bg=biograph(GTM,ids,'LayoutScale',0.3);
colorEdges(bg,ids,mat2);%bg.showWeights='On';
g=biograph.bggui(bg);
fig2=figure;a=gca;
copyobj(allchild(g.biograph.hgAxes),a);%Copy biograph viewer to current axes
set(0,'ShowHiddenHandles','on');close(setdiff(get(0,'Children'),[fig fig2]));%Close the biograph viewer
axes(a);
axis off
set(gcf,'Color','w');
end

%Return result
fulltree=combine(ttree,ptreeFromFullTree(record(end).tree));

function colorEdges(bg,ids,mat)
for i=1:size(mat,1)
    for j=1:size(mat,1)
        GroupA = ids(i);
        GroupB = ids(j);
        edgesSel = getedgesbynodeid(bg,[GroupA;GroupB]');
        set(edgesSel,'LineColor',[1 1 1]*max(0,min(1,(1-mat(i,j)))));
        set(edgesSel,'Label','lol');
    end
end
