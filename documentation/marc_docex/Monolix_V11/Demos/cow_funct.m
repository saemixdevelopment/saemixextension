function [f,g]=pat_funct(phi,x,id);
%
%   phi : individual parameters
%     x  : regression variables
%  
%%%%%%%%%%%%%%%%%%%%%%%
a=phi(id,1,:);
b=phi(id,2,:);
k=phi(id,3,:);
%%%%%%%%%%%%%%%%%%%%%%%

f=a.*(1-(b/1000).*exp(-k.*(x/100000)));
g=1;

