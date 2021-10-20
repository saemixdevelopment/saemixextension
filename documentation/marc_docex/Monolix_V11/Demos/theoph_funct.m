function  [f,g]=theoph_funct(phi,x,id);
%   phi : individual parameters
%     x  : regression variables
%%%%%%%%%%%%%%%%%%%%%%%
d=x(:,1,:);
t=x(:,2,:);
%%%%%%%%%%%%%%%%%%%%%%%
ke=exp(phi(id,1,:));
ka=exp(phi(id,2,:));
Cl=exp(phi(id,3,:));
%%%%%%%%%%%%%%%%%%%%%%%

f=d.*ka.*ke./(Cl.*(ka-ke)).*(exp(-ke.*t)-exp(-ka.*t));

g=1;
%g=d;
%g=f;
%g=(0.7 + f.^0.3);