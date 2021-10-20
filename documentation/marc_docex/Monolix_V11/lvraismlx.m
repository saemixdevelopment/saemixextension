function log_likelihood=lvraismlx(s);

%LVRAISMLX  computation of the log-likelihood, called by MONOLIX
%
%  S is a structure array that contains the information required for running
%  the algorithms.

%  Marc Lavielle
%  Version 1.1 ;  2005/02/18

%This program is free software; you can redistribute it and/or modify it
%under the terms of the GNU General Public License as published by the Free
%Software Foundation; either version 2 of the License, or (at your option)
%any later version.
%This program is distributed in the hope that it will be useful, but
%WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
%for more details (http://www.gnu.org/copyleft/gpl.html).


load(s.resultats)
K=s.Klog;
N=length(nb_measures);
nt_max=max(nb_measures);
nphi1=length(i1);
nphi0=length(i0);
nphi=nphi1+nphi0;

mprior_phi0=mprior_phi(:,i0);
mprior_phi1=mprior_phi(:,i1);
Gamma2_phi1=cov_random(i1,i1);

mtild_phi1=repmat(mpost_phi(:,i1),[1,1,K]);
vpost_phi1=vpost_phi(:,i1);
vtild_phi1=repmat(vpost_phi1,[1,1,K]);
%
phi1=mtild_phi1+sqrt(vtild_phi1).*randn(N,nphi1,K);
%
vprior_phi1=diag(Gamma2_phi1);
if max(max(abs(diag(vprior_phi1)-Gamma2_phi1)))>0
    IGamma2_phi1=inv(Gamma2_phi1);
    d1= phi1-repmat(mprior_phi1,[1,1,K]);
    d1=reshape(permute(d1,[2 1 3]),[nphi1,N*K,1])';
    d2=d1*IGamma2_phi1;
    d3=sum(d1.*d2,2);
    clear d1 d2
    e2=-reshape(d3,N,K)/2-0.5*(log(det(Gamma2_phi1)) + nphi1*log(2*pi)) ;
else
    pi_phi1=-sum(( (phi1-repmat(mprior_phi1,[1,1,K])).^2 )./repmat(vprior_phi1',[N,1,K]),2)/2;
    e2=reshape(pi_phi1,N,K)-0.5*(sum(log(vprior_phi1))+nphi1*log(2*pi));
end;
%
pitild_phi1=-sum(( (phi1-mtild_phi1).^2 )./vtild_phi1,2)/2;
e3=reshape(pitild_phi1,N,K)-(1/2).*repmat(sum(log(2*pi*vpost_phi1),2),1,K);
clear mtild_phi1 vtild_phi1 pitild_phi1 
%
phi=zeros(N,nphi,K);
phi(:,i1,:)=phi1;
phi(:,i0,:)=repmat(mpost_phi(:,i0),[1,1,K]);
clear phi1
t=repmat(design,[1,1,K]);
[f,g]=feval(fonction,phi,t,id);
f=squeeze(f);
g=squeeze(g);
Y=repmat(y,1,K);
DYF=zeros(N*nt_max,K);
DYF(indio,:)=((Y-f)./g).^2;
DYF=reshape(DYF,[nt_max,N,K]);
DG=zeros(N*nt_max,K);
DG(indio,:)=log(2*pi*sigma2*(g.^2))/2;
DG=reshape(DG,[nt_max,N,K]);

eyf=squeeze(sum(DYF,1));
eg=squeeze(sum(DG,1));

clear Y f  g DYF DG
e1=-eyf./(2*sigma2)-eg;
%
a=e1+e2-e3;
log_likelihood=sum(log(mean(exp(a),2)));
eval(['save ',s.resultats,' log_likelihood -append']);

