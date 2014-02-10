function plotTTree(ttree,method)
%Plots a transmission tree
if nargin==1,method=2;end
n=size(ttree,1);
a=gca;
fig=gcf;

if method==1
    %First plotting method, using one line per host
    hold on;
    for i=1:n
        plot([ttree(i,1) ttree(i,2)],[i i]);
        if ttree(i,3)>0,plot([ttree(i,1) ttree(i,1)],[ttree(i,3) i]);end
    end
    ylim([0 n+1]);
end

if method==2
    %Second plotting method, using a biograph
    %if n==1,return;end%The biograph does not work with only one node
    cm=zeros(n+1,n+1);
    ids{1}='SOURCE';
    for i=1:n
        cm(ttree(i,3)+1,i+1)=1;
        ids{i+1}=sprintf('%d [%.2f;%.2f]',i,ttree(i,1),ttree(i,2));
        %ids{i}=sprintf('%d',i);
end
    bg=biograph(cm,ids,'LayoutScale',0.5);%Create biograph
    set(bg.nodes(:),'FontSize',12);
    g=biograph.bggui(bg);%View biograph
    copyobj(allchild(g.biograph.hgAxes),a);%Copy biograph viewer to current axes
    set(0,'ShowHiddenHandles','on');close(setdiff(get(0,'Children'),fig));%Close the biograph viewer
    axes(a);
    axis off
    set(gcf,'Color','w');
end