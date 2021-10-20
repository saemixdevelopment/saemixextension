function [f,g]=yieldQP_funct(phi,x,id);
%
%   phi : individual parameters
%     x  : regression variables
%  
%%%%%%%%%%%%%%%%%%%%%%%
Y_max=phi(id,1,:);
X_max=phi(id,2,:);
B=phi(id,3,:);
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
f=Y_max;
i=find(x<X_max);
f(i)=Y_max(i)+B(i).*(x(i)-X_max(i)).^2;
g=1;
%%%%%%%%%%%%%%%%%%%%%%%

