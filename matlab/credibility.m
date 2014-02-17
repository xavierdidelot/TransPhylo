function r=credibility(vals,perc)
if nargin==1,perc=0.95;end
v=sort(vals);
l=length(v);
min=v(floor(l*(1-perc)/2));
max=v(floor(l*((1-perc)/2+perc)));
r=[mean(v) min max];