function [f,g]=PD_funct(phi,x,id);
%
%   phi : individual parameters
%     x  : regression variables
%  
%%%%%%%%%%%%%%%%%%%%%%%
E0=exp(phi(id,1,:));
Emax=exp(phi(id,2,:));
ED50=exp(phi(id,3,:));
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
DOSE=x(:,1,:);
%%%%%%%%%%%%%%%%%%%%%%%
f=E0+Emax.*DOSE./(ED50+DOSE);
g=1;
%%%%%%%%%%%%%%%%%%%%%%%


