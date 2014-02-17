function plotBothTree(tree,method,names)
%Plot both phylogenetic and transmission trees on the same figure
%The second parameter indicates whether to show stars for transmission
%events (method=0), color branches according to host (method=1) or use
%boxes to represent hosts (method=2)
n=sum(tree(:,2)+tree(:,3)==0);
if nargin<3,for i=1:n,names{i}=sprintf('%d',i);end;end
if nargin<2,method=1;end
hold on
host=tree(:,4);
cmap=colormap;

%Assign a unique color to each host
cols=zeros(n,3);
for i=1:n
    cols(i,:)=cmap(min(size(cmap,1),1+floor(size(cmap,1)*i/n)),:);
end

%Assign alternating colors to hosts
%basiccols=cmap(1:round(size(cmap,1)/6):size(cmap,1),:);
%basiccols=[1 0 0;0 1 0;0 0 1];
%co=zeros(size(tree,1),1);
%td=size(tree,1);
%co(td)=0;
%while ~isempty(td)
%    chi=setdiff(tree(td(1),2:3),0);
%    if length(chi)==1,co(chi)=1+mod(co(td(1)),size(basiccols,1));else co(chi)=co(td(1));end
%    td=[td(2:end) chi];
%end
%cols=basiccols(co(1:n),:);

%Determine ys
ys=zeros(n,1);
todo=[size(tree,1) 0 0.5 1];%Matrix of nodes to do, with associated starting x and y coordinates and scale
while size(todo,1)>0
    w=todo(1,1);x=todo(1,2);y=todo(1,3);scale=todo(1,4);
    if tree(w,2)==0 && tree(w,3)==0
        %Leaf node
        ys(w)=y;
    elseif tree(w,3)==0
        %Transmission node
        todo=[todo;tree(w,2) tree(w,1) y scale];
    else
        %Binary node
        todo=[todo;tree(w,2) tree(w,1) y+scale/2 scale/2;tree(w,3) tree(w,1) y-scale/2 scale/2];
    end
    todo=todo(2:end,:);%Remove first node
end
[~,ys]=sort(ys);
ys(ys)=1:length(ys);
ylims=[ys ys];
for i=n+1:size(tree,1)
    f=1+find(tree(i,2:3)>0);
    ylims(i,1)=min(ylims(tree(i,f),1));
    ylims(i,2)=max(ylims(tree(i,f),2));
    ys(i)=mean(ys(tree(i,f)));
end

todo=[size(tree,1) tree(end,1)];%Matrix of nodes to do, with associated starting x
while size(todo,1)>0
    w=todo(1,1);x=todo(1,2);y=ys(w);
    if method==1&&host(w)>0,col=cols(host(w),:);else
    if tree(w,2)>0,col=cols(host(tree(w,2)),:);else col=[0 0 0];end;end
    if tree(w,2)==0 && tree(w,3)==0
        %Leaf node
        if method==2
           patch([tree(w,1) tree(w,1) max(tree(:,1)) max(tree(:,1))],[ylims(w,1)-0.5 ylims(w,2)+0.5 ylims(w,2)+0.5 ylims(w,1)-0.5],[1 1 1],'LineStyle','none');
        end
        p=plot([x tree(w,1)],[y y],'Color',col*(method==1),'LineWidth',2);
        uistack(p,'bottom');
        text(tree(w,1)+(max(tree(:,1))-min(tree(:,1)))/100,y,names{w},'FontSize',12);
    elseif tree(w,3)==0
        %Transmission node
        p=plot([x tree(w,1)],[y y],'Color',col*(method==1),'LineWidth',2);
        uistack(p,'bottom');
        if method==0||method==1
            plot(tree(w,1),y,'r*','MarkerSize',10);
        elseif method==2
            patch([tree(w,1) tree(w,1) max(tree(:,1)) max(tree(:,1))],[ylims(w,1)-0.5 ylims(w,2)+0.5 ylims(w,2)+0.5 ylims(w,1)-0.5],col,'LineStyle','none');
%            patch([tree(fi,1) tree(fi,1) max(tree(:,1)) max(tree(:,1))],[y-scale y+scale y+scale y-scale],[1 1 1],'LineStyle','none');
        end
        todo=[todo;tree(w,2) tree(w,1)];
    else
        %Binary node
        p=plot([x tree(w,1)],[y y],'Color',col*(method==1),'LineWidth',2);
        uistack(p,'bottom');
        p=plot([tree(w,1) tree(w,1)],[ys(tree(w,2)) ys(tree(w,3))],'Color',col*(method==1),'LineWidth',2);
        uistack(p,'bottom');
        todo=[todo;tree(w,2) tree(w,1);tree(w,3) tree(w,1)];
    end
    todo=todo(2:end,:);%Remove first node
end
set(gcf,'Color','w');
set(gca,'YTick',[]);
set(gca,'YColor','w');
ylim([0 n+1]);