clear all
close all
rand('state',0);
neg=5;
for i=1:10000
    times=10*ones(10,1);
    nodes=withinhost(times,neg);
    tmrca(i)=max(nodes(1:end-1,1))-min(nodes(1:end-1,1));
    tmrca2(i)=11;
    while tmrca2(i)>10
        tmrca2(i)=0;
        for j=2:length(times),tmrca2(i)=tmrca2(i)-log(rand)/nchoosek(j,2);end
        tmrca2(i)=tmrca2(i)*neg;
    end
end
subplot(2,1,1);
title('Distribution of TMRCA using withinhost code');
hist(tmrca,0.5:10.5);
subplot(2,1,2);
title('Distribution of TMRCA using standard coalescent');
hist(tmrca2,0.5:10.5);