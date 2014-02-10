function sane=sanitycheck(fulltree)
ttree=ttreeFromFullTree(fulltree);
sane=1;
for i=1:size(ttree,1)
    infector=ttree(i,3);
    if infector==0,continue;end
    if ttree(i,1)>ttree(infector,2),sprintf('Transmission after recovery of infector for %d->%d, %f>%f',infector,i,ttree(i,1),ttree(infector,2)),sane=0;end
    if ttree(i,1)<ttree(infector,1),sprintf('Transmission before infection of infector for %d->%d, %f<%f',infector,i,ttree(i,1),ttree(infector,2)),sane=0;end
end