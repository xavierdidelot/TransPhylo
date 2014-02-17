function plotPTree2(tree,names)
%Takes a full tree as input
n=sum(tree(:,2)+tree(:,3)==0);
if nargin<2,for i=1:n,names{i}=sprintf('%d',i);end;end
hold on
host=tree(:,4);

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

todo=[size(tree,1)-1 tree(end-1,1)];%Matrix of nodes to do, with associated starting x
while size(todo,1)>0
    w=todo(1,1);x=todo(1,2);y=ys(w);
    col=[0 0 0];
    %if method==1&&host(w)>0,col=cols(host(w),:);else
    %if tree(w,2)>0,col=cols(host(tree(w,2)),:);else col=[0 0 0];end;end
    if tree(w,2)==0 && tree(w,3)==0
        %Leaf node
        %if method==2
        %   patch([tree(w,1) tree(w,1) max(tree(:,1)) max(tree(:,1))],[ylims(w,1)-0.5 ylims(w,2)+0.5 ylims(w,2)+0.5 ylims(w,1)-0.5],[1 1 1],'LineStyle','none');
        %end
        p=plot([x tree(w,1)],[y y],'Color',col,'LineWidth',2);
        uistack(p,'bottom');
        text(tree(w,1)+(max(tree(:,1))-min(tree(:,1)))/100,y,names{w},'FontSize',12);
    elseif tree(w,3)==0
        %Transmission node
        p=plot([x tree(w,1)],[y y],'Color',col,'LineWidth',2);
        uistack(p,'bottom');
        %if method==0||method==1
        %    plot(tree(w,1),y,'r*','MarkerSize',10);
        %elseif method==2
        %    patch([tree(w,1) tree(w,1) max(tree(:,1)) max(tree(:,1))],[ylims(w,1)-0.5 ylims(w,2)+0.5 ylims(w,2)+0.5 ylims(w,1)-0.5],col,'LineStyle','none');
%            patch([tree(fi,1) tree(fi,1) max(tree(:,1)) max(tree(:,1))],[y-scale y+scale y+scale y-scale],[1 1 1],'LineStyle','none');
        %end
        todo=[todo;tree(w,2) tree(w,1)];
    else
        %Binary node
        p=plot([x tree(w,1)],[y y],'Color',col,'LineWidth',2);
        uistack(p,'bottom');
        p=plot([tree(w,1) tree(w,1)],[ys(tree(w,2)) ys(tree(w,3))],'Color',col,'LineWidth',2);
        uistack(p,'bottom');
        todo=[todo;tree(w,2) tree(w,1);tree(w,3) tree(w,1)];
    end
    todo=todo(2:end,:);%Remove first node
end
set(gcf,'Color','w');
set(gca,'YTick',[]);
set(gca,'YColor','w');
ylim([0 n+1]);