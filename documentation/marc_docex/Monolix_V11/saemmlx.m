function algosaem(s,hfig)
%SAEMMLX  SAEM algorithm, called by MONOLIX
%
%  S is a structure array that contains the information required for running
%  the algorithms.
%
%  ALGOSAEM(S) compute the maximum likelihood estimator of the parameters
%  of the model, the standard deviations and the p-value obtained from the
%  Wald test.
%
%  ALGOSAEM(S,HFIG)  plots the estimations in the figure handle HFIG.

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

tic
fonction=getfield(s,'fonction');
donnees=getfield(s,'donnees');
resultats=getfield(s,'resultats');
mcov=getfield(s,'mcov');
lcov_ini=getfield(s,'lcov_ini');
var_phi_ini=getfield(s,'var_phi_ini');
va=getfield(s,'va');
Klog=getfield(s,'Klog');
minv=getfield(s,'minv');
kitaff=getfield(s,'kitaff');
nu=getfield(s,'nu');
rmcmc=getfield(s,'rmcmc');
nmc=getfield(s,'nmc');
vna=getfield(s,'vna');
niter0=getfield(s,'niter0');
sigma2_ini=getfield(s,'sigma2_ini');
covstruct=getfield(s,'covstruct');
nb_varex=getfield(s,'nb_varex');
nb_param=getfield(s,'nb_param');
nb_covariates=getfield(s,'nb_covariates');
logstruct=getfield(s,'logstruct');
iop_var=getfield(s,'iop_var');
coef_sa=getfield(s,'coef_sa');
nb_sa     =vna(1)*.75;
nb_correl=vna(1)*.75;
minv=max(minv,1e-20);
%
[ncov,nphi]=size(mcov);
Mstr={};
Mstra={};
Mstrb={};
l=0;
for j=1:nphi
    if mcov(1,j)==1;
        l=l+1;
        Mstr{l}=['e^{\mu_' num2str(j),'}'];
        Mstra{l}=['mu_' num2str(j)];
        Mstrb{l}=['exp(mu_' num2str(j),')'];
    end;
    if ncov==2 & mcov(2,j)==1;
        l=l+1;
        Mstr{l}=['\beta_' num2str(j)];
        Mstrb{l}=['beta_' num2str(j)];
    end;
    if ncov>2
        for k=2:ncov
            if mcov(k,j)==1;
                l=l+1;
                Mstr{l}=['\beta_' num2str(j) '_' num2str(k-1)];
                Mstrb{l}=['beta_' num2str(j) num2str(k-1)];
            end;
        end;
    end;
end;
for j=1:nphi
    if covstruct(j,j)==1;
        l=l+1;
        if iop_var==1
            Mstr{l}=['\omega^2_' num2str(j)];
            Mstrb{l}=['omega2_' num2str(j)];
        else
            Mstr{l}=['\omega_' num2str(j)];
            Mstrb{l}=['omega_' num2str(j)];
        end;
    end;
end;
l=l+1;
if iop_var==1
    Mstr{l}='\sigma^2';
    Mstrb{l}='sigma2';
else
    Mstr{l}='\sigma';
    Mstrb{l}='sigma';
end;
%P=[10,50,800,500];
P=[.5 .44 .5 .5];
if nargin>1
    figure(hfig);
    clf
    addmenumlx(hfig);
    set(hfig,'name',' Estimation of the parameters');
    set(hfig,'unit','normalized','position',P);
    set(hfig,'visible','on');
end;
format short g

%eval(['load ' param_data])
%%%
[y,covariables,design,io,id,nb_measures]=datasaem(donnees,nb_varex);
[N,nt_max]=size(io);
[N,nb_covariables]=size(covariables);
ntotal=sum(nb_measures);
%
yM=repmat(y,nmc,1);
designM=repmat(design,nmc,1);
ioM=repmat(io,nmc,1);
idM=repmat(id,nmc,1)+kron(N*(0:nmc-1)',ones(ntotal,1));
%nb_measuresM=repmat(nb_measures,nmc,1);
%covariablesM=repmat(covariables,nmc,1);
nM=N*nmc;
%
indio=find(io');
DYF=zeros(size(ioM))';
indioM=find(ioM');
niter=sum(vna);
niter_phi0=round(niter/2);
coef_phi0=0.0001^(1/niter_phi0);
coefvar_phi0_ini=0.1;
pash=[ones(1,niter0) 1./(1:niter)];
%covy=cov(y(:));
%sigma2_ini=var_sigma2_ini;
niter=sum(vna);

Mcovariables=[ones(N,1) covariables];
%
jlog1=find(logstruct==1);
jcov=find(sum(mcov,2)>0);
Mcovariables=Mcovariables(:,jcov);
mcov=mcov(jcov,:);
lcov_ini=lcov_ini(jcov,:);
lcov_ini(1,jlog1)=log(lcov_ini(1,jlog1));
j0=find(mcov==0);
lcov_ini(j0)=0;
%
%DGamma2_phi_ini=max((mean(Mcovariables*lcov_ini,1).*var_phi_ini).^2,1);
DGamma2_phi_ini=var_phi_ini;
i1=find(diag(covstruct)~=0);
i0=find(diag(covstruct)==0);
covstruct1=covstruct(i1,i1);
DGamma2_phi1_ini=DGamma2_phi_ini(i1);
Gamma2_phi1_ini=diag(DGamma2_phi1_ini);
DGamma2_phi0_ini=(mean(Mcovariables*lcov_ini(:,i0),1)*coefvar_phi0_ini).^2;
Gamma2_phi0_ini=diag(DGamma2_phi0_ini);
nphi1=length(i1);
nphi0=length(i0);
%
[ncov,nphi]=size(mcov);
mcov1=mcov(:,i1);
mcov0=mcov(:,i0);
ind_cov=find(mcov==1);
ind_cov1=find(mcov1==1);
ind_cov0=find(mcov0==1);
lcov1_ini=lcov_ini(:,i1);
lcov0_ini=lcov_ini(:,i0);
%
na=length(va);
pas=1./(1:vna(1)).^va(1) ;
for ia=2:na
    k1=pas(end)^(-1/va(ia));
    pas=[pas 1./(k1+1:k1+vna(ia)).^va(ia) ];
end;
%%%%%%
nlambda1=sum(sum(mcov1==1));
nlambda0=sum(sum(mcov0==1));
nlambda=nlambda1+nlambda0;
nd1=nphi1+nlambda1+1;
nd2=nphi1+nlambda1+nlambda0;
nb_param=nd2+1;
%
Plambda_ini=lcov_ini(ind_cov);
Plambda_ini=Plambda_ini(:)';
%
Plambda1_ini=zeros(1,0);
mprior_phi1=zeros(N,nphi1);
pc=sum(mcov,1);
ipc= cumsum([0 pc(1:nphi-1)])+1;
ipcl1=ipc(jlog1);
pc1=zeros(1,nphi1);
COV1=zeros(N,0);
LCOV1=zeros(nlambda1,nphi1);
MCOV1=zeros(nlambda1,nphi1);
j1=1;
for j=1:nphi1
    j_cov=find(mcov1(:,j)==1);
    lambdaj=lcov1_ini(j_cov,j);
    Aj=Mcovariables(:,j_cov);
    COV1=[COV1 Aj];
    nlj=length(lambdaj);
    j2=j1+nlj-1;
    LCOV1(j1:j2,j)=ones(nlj,1);
    j1=j2+1;
    mprior_phi1(:,j)=Aj*lambdaj;
    Plambda1_ini=[Plambda1_ini lambdaj'];
    pc1(j)=length(lambdaj);
end;
COV21=COV1'*COV1;
jcov1=find(LCOV1==1);
MCOV1(jcov1)=Plambda1_ini;
%
Plambda0_ini=zeros(1,0);
mprior_phi0=zeros(N,nphi0);
pc0=zeros(1,nphi0);
COV0=zeros(N,0);
LCOV0=zeros(nlambda0,nphi0);
MCOV0=zeros(nlambda0,nphi0);
j1=1;
for j=1:nphi0
    j_cov=find(mcov0(:,j)==1);
    lambdaj=lcov0_ini(j_cov,j);
    Aj=Mcovariables(:,j_cov);
    COV0=[COV0 Aj];
    nlj=length(lambdaj);
    j2=j1+nlj-1;
    LCOV0(j1:j2,j)=ones(nlj,1);
    j1=j2+1;
    mprior_phi0(:,j)=Aj*lambdaj;
    Plambda0_ini=[Plambda0_ini lambdaj'];
    pc0(j)=length(lambdaj);
end;
COV20=COV0'*COV0;
jcov0=find(LCOV0==1);
MCOV0(jcov0)=Plambda0_ini;
%
indiphi=[];
for k=1:nphi1
    indiphi=[indiphi repmat(i1(k),1,pc1(k))];
end;
for k=1:nphi0
    indiphi=[indiphi repmat(i0(k),1,pc0(k))];
end;
[temp,indiphi]=sort(indiphi);
%
PARA=zeros(niter+1,nb_param);
Plambda_ini(ipcl1)=exp(Plambda_ini(ipcl1));
if iop_var==1
    PARA(1,:)= [ Plambda_ini diag(Gamma2_phi1_ini)'  sigma2_ini ];
else
    PARA(1,:)= [ Plambda_ini sqrt(diag(Gamma2_phi1_ini)')  sqrt(sigma2_ini) ];
end;
%   initialisation des parametres
sigma2=sigma2_ini;
Gamma2_phi1=Gamma2_phi1_ini;
Gamma2_phi0=Gamma2_phi0_ini;
Gamma_phi1=sqrtm(Gamma2_phi1);
Gamma_phi0=sqrtm(Gamma2_phi0);
phi=zeros(N,nphi,nmc);
for k=1:nmc
    phi(:,i1,k)=randn(N,nphi1)*Gamma_phi1+mprior_phi1;
    phi(:,i0,k)=randn(N,nphi0)*Gamma_phi0+mprior_phi0;
end;
phiM=zeros(nM,nphi);
phiM(:,i1)=randn(nM,nphi1)*Gamma_phi1+repmat(mprior_phi1,nmc,1);
phiM(:,i0)=randn(nM,nphi0)*Gamma_phi0+repmat(mprior_phi0,nmc,1);
%   initialisation des statistiques
statphi11=sum(phi(:,i1,:),3)/nmc;
statphi12=0;
statphi01=sum(phi(:,i0,:),3)/nmc;
statphi02=0;
staty=0;
statfg1=0;
statfg2=0;
resy=zeros(nmc,1);
for k=1:nmc
    phik=phi(:,:,k);
    phi1k=phik(:,i1);
    phi0k=phik(:,i0);
    statphi12=statphi12+phi1k'*phi1k;
    statphi02=statphi02+phi0k'*phi0k;
    [fk,gk]=feval(fonction,phik,design,id);
    resy(k)=sum(((y-fk)./gk).^2);
end;
statphi12=statphi12/nmc;
statphi02=statphi02/nmc;
statresy=sum(resy)/nmc;
Delta=0;
G=0;
Ha=0;
Hb=0;
L=0;
mpost_phi=0;
cpost_phi=0;
LOGL=[];
%%%%%%%%%%%%%%%%  ALGO  %%%%%%%%%%%%%%%%%%
for kiter=1:niter;
    gamma2_phi1=diag(Gamma2_phi1);
    IGamma2_phi1=inv(Gamma2_phi1);
    D1Gamma21=LCOV1*IGamma2_phi1;
    D2Gamma21=D1Gamma21*LCOV1';
    CGamma21=COV21.*D2Gamma21;
    gamma2_phi0=diag(Gamma2_phi0);
    IGamma2_phi0=inv(Gamma2_phi0);
    D1Gamma20=LCOV0*IGamma2_phi0;
    D2Gamma20=D1Gamma20*LCOV0';
    CGamma20=COV20.*D2Gamma20;
    %    algo  MCMC
    Gamma_phi1=sqrtm(Gamma2_phi1);
    Gdiag_phi1=sqrtm(diag(diag(Gamma2_phi1)))*rmcmc;
    Gamma_phi0=sqrtm(Gamma2_phi0);
    Gdiag_phi0=sqrtm(diag(diag(Gamma2_phi0)))*rmcmc;
    %
    %  MCMC
    %
    if kiter==1
        nu1=20*nu(1);
        nu2=20*nu(2);
        nu3=20*nu(3);
    else
        nu1=nu(1);
        nu2=nu(2);
        nu3=nu(3);
    end;
    [f,g]=feval(fonction,phiM,designM,idM);
    DYF(indioM)=0.5/sigma2*((yM-f)./g).^2+log(g);
    U_y=sum(DYF,1)';
    phiMc=phiM;
    mprior_phi1M=repmat(mprior_phi1,nmc,1);
    mprior_phi0M=repmat(mprior_phi0,nmc,1);
    if nphi1>0
        for u=1:nu1
            phiMc(:,i1)=randn(nM,nphi1)*Gamma_phi1+mprior_phi1M;
            [fc,gc]=feval(fonction,phiMc,designM,idM);
            DYF(indioM)=0.5/sigma2*((yM-fc)./gc).^2+log(gc);
            Uc_y=sum(DYF,1)';
            deltu=Uc_y-U_y;
            ind=find( deltu<-log(rand(nM,1)) );
            phiM(ind,i1)=phiMc(ind,i1);
            U_y(ind,:)=Uc_y(ind,:);
        end;
        dphi=phiM(:,i1)-mprior_phi1M;
        U_phi=0.5*sum(dphi.*(dphi*IGamma2_phi1),2);
        for u=1:nu2
            phiMc(:,i1)=phiM(:,i1)+randn(nM,nphi1)*Gdiag_phi1;
            [fc,gc]=feval(fonction,phiMc,designM,idM);
            DYF(indioM)=0.5/sigma2*((yM-fc)./gc).^2+log(gc);
            Uc_y=sum(DYF,1)';
            dphic=phiMc(:,i1)-mprior_phi1M;
            Uc_phi=0.5*sum(dphic.*(dphic*IGamma2_phi1),2);
            deltu=Uc_y-U_y+Uc_phi-U_phi;
            ind=find( deltu<-log(rand(nM,1)) );
            phiM(ind,i1)=phiMc(ind,i1);
            U_y(ind,:)=Uc_y(ind,:);
            U_phi(ind,:)=Uc_phi(ind,:);
        end;
        for u=1:nu3
            for k1=1:nphi1
                phiMc=phiM;
                phiMc(:,i1(k1))=phiM(:,i1(k1))+randn(nM,1)*Gdiag_phi1(k1,k1);
                [fc,gc]=feval(fonction,phiMc,designM,idM);
                DYF(indioM)=0.5/sigma2*((yM-fc)./gc).^2+log(gc);
                Uc_y=sum(DYF,1)';
                dphic=phiMc(:,i1)-mprior_phi1M;
                Uc_phi=0.5*sum(dphic.*(dphic*IGamma2_phi1),2);
                deltu=Uc_y-U_y+Uc_phi-U_phi;
                ind=find( deltu<-log(rand(nM,1)) );
                phiM(ind,i1)=phiMc(ind,i1);
                U_y(ind,:)=Uc_y(ind,:);
                U_phi(ind,:)=Uc_phi(ind,:);
            end;
        end;
    end;
    %%%   simulation de phi0
    if nphi0>0
        phiMc=phiM;
        t=cputime;
        for u=1:nu1
            phiMc(:,i0)=randn(nM,nphi0)*Gamma_phi0+mprior_phi0M;
            [fc,gc]=feval(fonction,phiMc,designM,idM);
            DYF(indioM)=0.5/sigma2*((yM-fc)./gc).^2+log(gc);
            Uc_y=sum(DYF,1)';
            deltu=Uc_y-U_y;
            ind=find( deltu<-log(rand(nM,1)) );
            phiM(ind,i0)=phiMc(ind,i0);
            U_y(ind,:)=Uc_y(ind,:);
        end;
        dphi=phiM(:,i0)-mprior_phi0M;
        U_phi=0.5*sum(dphi.*(dphi*IGamma2_phi0),2);
        for u=1:nu2
            phiMc(:,i0)=phiM(:,i0)+randn(nM,nphi0)*Gdiag_phi0;
            [fc,gc]=feval(fonction,phiMc,designM,idM);
            DYF(indioM)=0.5/sigma2*((yM-fc)./gc).^2+log(gc);
            Uc_y=sum(DYF,1)';
            dphic=phiMc(:,i0)-mprior_phi0M;
            Uc_phi=0.5*sum(dphic.*(dphic*IGamma2_phi0),2);
            deltu=Uc_y-U_y+Uc_phi-U_phi;
            ind=find( deltu<-log(rand(nM,1)) );
            phiM(ind,i0)=phiMc(ind,i0);
            U_y(ind,:)=Uc_y(ind,:);
            U_phi(ind,:)=Uc_phi(ind,:);
        end;
        for u=1:nu3
            for k0=1:nphi0
                phiMc=phiM;
                phiMc(:,i0(k0))=phiM(:,i0(k0))+randn(nM,1)*Gdiag_phi0(k0,k0);
                [fc,gc]=feval(fonction,phiMc,designM,idM);
                DYF(indioM)=0.5/sigma2*((yM-fc)./gc).^2+log(gc);
                Uc_y=sum(DYF,1)';
                dphic=phiMc(:,i0)-mprior_phi0M;
                Uc_phi=0.5*sum(dphic.*(dphic*IGamma2_phi0),2);
                deltu=Uc_y-U_y+Uc_phi-U_phi;
                ind=find( deltu<-log(rand(nM,1)) );
                phiM(ind,i0)=phiMc(ind,i0);
                U_y(ind,:)=Uc_y(ind,:);
                U_phi(ind,:)=Uc_phi(ind,:);
            end;
        end;
    end;
    %
    for k=1:nmc
        phi(:,:,k)=phiM((k-1)*N+1:k*N,:);
    end;
    %    approximation stochastique
    Statphi11=sum(phi(:,i1,:),3);
    Statphi12=0;
    Statphi01=sum(phi(:,i0,:),3);
    Statphi02=0;
    Staty=0;
    Statfg1=0;
    Statfg2=0;
    resy=zeros(nmc,1);
    D1=0;
    D11=0;
    D2=0;
    d2logk=zeros(nb_param,nb_param);
    d2logk(1:nlambda1,1:nlambda1)=-CGamma21;
    d2logk(nlambda1+1:nlambda,nlambda1+1:nlambda)=-CGamma20;
    %d2logk(nd1:nd2,nd1:nd2)=-IGamma2_phi0;
    for k=1:nmc
        phik=phi(:,:,k);
        phi1k=phik(:,i1);
        phi0k=phik(:,i0);
        Statphi12=Statphi12+phi1k'*phi1k;
        Statphi02=Statphi02+phi0k'*phi0k;
        [fk,gk]=feval(fonction,phik,design,id);
        Staty=Staty+(y./gk)'*(y./gk);
        Statfg1=Statfg1+(y./gk)'*(fk./gk);
        Statfg2=Statfg2+(fk./gk)'*(fk./gk);
        resy(k)=sum(((y-fk)./gk).^2);
        %
        dphi1k=phi1k-mprior_phi1;
        dphi0k=phi0k-mprior_phi0;
        sdg1=sum(dphi1k.^2,1)'./gamma2_phi1;
        Md1=(IGamma2_phi1*(dphi1k'*Mcovariables))';
        Md0=(IGamma2_phi0*(dphi0k'*Mcovariables))';
        d1_mu_phi1=Md1(ind_cov1);
        d1_mu_phi0=Md0(ind_cov0);
        d1_loggamma2_phi1=0.5*sdg1-0.5*N;
        d1_logsigma2=0.5*resy(k)/sigma2-0.5*ntotal;
        d1logk=[d1_mu_phi1(:);d1_mu_phi0(:);d1_loggamma2_phi1;d1_logsigma2];
        D1=D1+d1logk;
        %Mdl1(kiter,:,k)=d1logk';
        D11=D11+d1logk*d1logk';
        %
        % v2phi=-sum(dphik,1)./gamma2_phi';

        w2phi=-0.5*sdg1';
        l=0;
        for j=1:nphi1
            for jj=1:pc1(j)
                l=l+1;
                temp=-COV1(:,l)'*dphi1k(:,j)/gamma2_phi1(j);
                d2logk(l,nlambda+j)=temp;
                d2logk(nlambda+j,l)=temp;
            end;
            d2logk(nlambda+j,nlambda+j)=w2phi(j);
        end;
        d2logk(nb_param,nb_param)=-0.5*resy(k)/sigma2;
        D2=D2+d2logk;
    end;
    %
    statphi11=statphi11+pas(kiter)*(Statphi11/nmc-statphi11);
    statphi12=statphi12+pas(kiter)*(Statphi12/nmc-statphi12);
    statphi01=statphi01+pas(kiter)*(Statphi01/nmc-statphi01);
    statphi02=statphi02+pas(kiter)*(Statphi02/nmc-statphi02);
    staty=staty+pas(kiter)*(Staty/nmc-staty);
    statfg1=statfg1+pas(kiter)*(Statfg1/nmc-statfg1);
    statfg2=statfg2+pas(kiter)*(Statfg2/nmc-statfg2);
    %    reactualisation des parametres
    Plambda1=(inv(CGamma21)*sum((D1Gamma21.*(COV1'*statphi11)),2))';
    Plambda0=(inv(CGamma20)*sum((D1Gamma20.*(COV0'*statphi01)),2))';
    MCOV1(jcov1)=Plambda1;
    MCOV0(jcov0)=Plambda0;
    mprior_phi1=COV1*MCOV1;
    mprior_phi0=COV0*MCOV0;
    %
    G1=statphi12/N+mprior_phi1'*mprior_phi1/N - statphi11'*mprior_phi1/N - mprior_phi1'*statphi11/N;
    %    Gamma2_phi1sa=Gamma2_phi1sa*coef_sa;
    %    Gamma2_phi1=max(Gamma2_phi1sa,G1)
    if kiter<=nb_sa
        Gamma2_phi1=max(Gamma2_phi1*(coef_sa),diag(diag(G1)));
    else
        Gamma2_phi1=G1;
    end;
    Gamma2_phi1=Gamma2_phi1.*covstruct1;
    Gmin=minv(i1);
    jDmin=find(diag(Gamma2_phi1)<Gmin');
    Gamma2_phi1(jDmin,jDmin)=Gmin(jDmin);
    if kiter<=nb_correl
        Gamma2_phi1=diag(diag(Gamma2_phi1));
    end;
    %
    if nphi0>0
        Gamma2_phi0=statphi02/N+mprior_phi0'*mprior_phi0/N - statphi01'*mprior_phi0/N - mprior_phi0'*statphi01/N;
        Gmin=minv(i0);
        jDmin=find(diag(Gamma2_phi0)<Gmin');
        Gamma2_phi0(jDmin,jDmin)=Gmin(jDmin);
        if kiter>niter_phi0
            dGamma2_phi0=dGamma2_phi0*coef_phi0;
        else
            dGamma2_phi0=diag(Gamma2_phi0);
        end;
        Gamma2_phi0=diag(dGamma2_phi0);
        %if Gamma2_phi0==NaN; keyboard; end;
        %disp(Gamma2_phi0)
    end;
    %
    statresy=staty-2*statfg1+statfg2;
    sig2=statresy/(ntotal);
    %    sig2sa=sig2sa*coef_sa;
    %    sigma2=max(sig2sa,sig2);
    if kiter<=nb_sa
        %      sigma2=max(sigma2*coef_sa,sig2);
        sigma2=sig2;
    else
        sigma2=sig2;
    end;
    %    information de Fisher
    %          Delta=Delta+pas(kiter)*(D1/nmc-Delta);
    %          G=G+pas(kiter)*(D11/nmc+D2/nmc-G);
    %          F_info=Delta*Delta'-G;
    L=L+pash(kiter)*(D1/nmc-L);
    DDa=(D1/nmc)*(D1/nmc)'-D11/nmc-D2/nmc;
    DDb=-D11/nmc-D2/nmc;
    Ha=Ha+pash(kiter)*(DDa- Ha);
    Hb=Hb+pash(kiter)*(DDb- Hb);
    mpost_phi=mpost_phi+pash(kiter)*(mean(phi,3)-mpost_phi);
    cpost_phi=cpost_phi+pash(kiter)*(mean(phi.^2,3)-cpost_phi);
    mpost_phi(:,i0)=mprior_phi0;
    %cpost_phi(:,i0)=0;
    %
    Plambda=[Plambda1 Plambda0];
    Plambda=Plambda(indiphi);
    ipc= cumsum([0 pc(1:nphi-1)])+1;
    PPP=Plambda;
    PPP(ipcl1)=exp(PPP(ipcl1));
    if iop_var==1
        PARA(kiter+1,:)= [ PPP diag(Gamma2_phi1)'  sigma2 ];
    else
        PARA(kiter+1,:)= [ PPP sqrt(diag(Gamma2_phi1)' ) sqrt(sigma2) ];
    end;
    if nargin>1 & kitaff>0
        if rem(kiter,kitaff)==0
            grafsaem(PARA,kiter,hfig,Mstr);
        end;
    end;
end;
clear PARA
F_info=Ha;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hphi=F_info(1:nlambda,:);
F_info(1:nlambda,:)=Hphi(indiphi,:);
Hphi=F_info(:,1:nlambda);
F_info(:,1:nlambda)=Hphi(:,indiphi);
D=sqrt(diag(inv(F_info)))';
%
vpost_phi=cpost_phi-mpost_phi.^2;
vpost_phi(:,i0)=0;
%
mprior_phi=zeros(N,nphi);
mprior_phi(:,i0)=mprior_phi0;
mprior_phi(:,i1)=mprior_phi1;
eval(['save ',resultats,' fonction y N nb_measures covariables design id  io indio']);
eval(['save ',resultats,' mprior_phi  sigma2  mpost_phi vpost_phi i1 i0 -append']);
eval(['save ',resultats,' F_info covstruct logstruct indiphi Mstr -append']);
%eval(['save ',resultats,' Plambda Gamma2_phi1 D -append']);
%%%%
P1=Plambda';
P2=log(diag(Gamma2_phi1));
P3=log(sigma2);
D1=D(1:nlambda)';
D2=D(nlambda+1:nlambda+nphi1)';
D3=D(nlambda+nphi1+1)';
%eP2=exp(P2+(D2.^2)/2);
eP2=exp(P2);
eD2=exp(P2).*sqrt(exp(2*(D2.^2))-exp(D2.^2));
%eP3=exp(P3+(D3.^2)/2);
eP3=exp(P3);
eD3=exp(P3).*sqrt(exp(2*(D3.^2))-exp(D3.^2));
%
Mstrm=str2mat(Mstrb);
j= find(Mstrm=='\');
Mstrm(j)=' ';
%
Mstrm1=Mstrm(1:nlambda,:);
Mstrm2=Mstrm(nlambda+1:nlambda+nphi1,:);
Mstrm3=Mstrm(nlambda+nphi1+1,:);
%
% res=s.resultats;
% j=find(res=='.');
% restxt=[res(1:j),'txt'];
% fid=fopen(restxt,'w');
%
disp([repmat(' ',1,size(Mstrm1,2)), '        estim. param.    stand. dev.   p-value'])
disp(' ')
if isreal(D1)
    pv1=2*(1-normcdf(abs(P1)./D1));
else
    D1=NaN*ones(size(D1));;
    pv1=NaN*ones(size(D1));;
end;
% P1(ipcl1)=exp(P1(ipcl1));
% D1(ipcl1)=P1(ipcl1).*sqrt(exp(2*(D1(ipcl1).^2))-exp(D1(ipcl1).^2));
% for k=1:nlambda
%     fprintf(fid,'%12s :  %12.3f %12.3f %12.3f \n',Mstrm1(k,:) , P1(k), D1(k) ,pv1(k));
%     disp(sprintf('%12s :  %12.3f %12.3f %12.3f ',Mstrm1(k,:) , P1(k), D1(k) ,pv1(k)));
% end;
R1=[P1, D1 ,pv1];
EP1=exp(P1(ipcl1));
ED1=EP1.*sqrt(exp(2*(D1(ipcl1).^2))-exp(D1(ipcl1).^2));
itestpc=zeros(1,nlambda);
itestpc(ipcl1)=1;
j=0;
for k=1:nlambda
    if itestpc(k)==1
        %        S1=[' log(',Mstr{k},')'];
        S1=Mstra{k};
    else
        S1=Mstrb{k};
        %                 S1=Mstrm1(k,:);
    end;
    %    fprintf(fid,'%12s :  %12.3g %12.3g %12.3g \n',S1, R1(k,:));
    disp(sprintf('%12s :  %12.3g %12.3g %12.3g ',S1 , R1(k,:)));
    if itestpc(k)==1
        j=j+1;
        %  fprintf(fid,'%12s :  %12.3g %12.3g \n',Mstrm1(ipcl1(j),:), EP1(j) , ED1(j));
        disp(sprintf('%12s :  %12.3g %12.3g  ',Mstrm(ipcl1(j),:) , EP1(j) , ED1(j)));
    end;
end;
P1(ipcl1)=EP1;
D1(ipcl1)=ED1;
disp(' ')
if isreal(eD2)
    pv2=(1-normcdf(abs(eP2)./eD2));
else
    D2=NaN*ones(size(eD2));;
    pv2=NaN*ones(size(eD2));;
end;
if iop_var==1
    R2=[eP2, eD2 ,pv2];
    %disp('  variance of the random effects:')
else
    R2=[sqrt(eP2), sqrt(exp(P2+D2.^2/2)-exp(P2+D2.^2/4)) ,pv2];
    %disp('  standard deviation of the random effects:')
end;
for k=1:nphi1
    disp(sprintf('%12s :  %12.3g %12.3g %12.3g',Mstrm2(k,:) , R2(k,:)));
    %   fprintf(fid,'%12s :  %12.3g %12.3g %12.3g \n',Mstrm2(k,:) , R2(k,:));
end;

covstruct1=covstruct(i1,i1);
testG=min(min(covstruct1==eye(nphi1)));
if testG==0
    d1=diag(Gamma2_phi1);
    sd1=sqrt(d1*d1');
    disp(' ')
    disp('  correlation matrix of the random effects:')
    format bank
    disp(Gamma2_phi1./sd1);
    format short g
    %   fprintf('  smallest eigenvalue of the correlation matrix: %0.3g\n',min(abs(eig(Gamma2_phi1./sd1))))
end;
disp(' ')
%disp('  residual variance:')
if isreal(eD3)==0
    eD3=NaN;
end;
if iop_var==1
    R3=[eP3, eD3 ];
    % disp('  residual variance:')
else
    R3=[sqrt(eP3) , sqrt(exp(P3+D3/2)-exp(P3+D3/4))];
    % disp('  residual standard deviation:')
end;
disp(sprintf('%12s :  %12.3g %12.3g  ',Mstrm3 , R3));
%fprintf(fid,'%12s :  %12.3g %12.3g  \n',Mstrm3 , R3);
%
%fclose(fid);
fixed_effects=P1;
sd_fixed=D1;

cov_random=zeros(nphi,nphi);
cov_random(i1,i1)=Gamma2_phi1;
pv_fixed=pv1;
sd_random=eD2;
pv_random=pv2;
sd_sigma2=eD3;
%eval(['save ',resultats,' eP2 eD2 eP3 eD3 -append']);
eval(['save ',resultats,' fixed_effects sd_fixed pv_fixed cov_random sd_random  pv_random sd_sigma2 -append']);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,covariables,design,indicobs,id,nb_measures]=datasaem(donnees,nb_varex);

S=load(donnees);
[Ss,ords]=sortrows(S,1);
[temp,ordr]=sort(ords);
ID=S(:,1);
[idu,itemp,ids]=unique(ID);
[temp,jtemp]=sort(itemp);
N=length(idu);
nb_covariables=size(S,2)-nb_varex-2;
nb_measures=zeros(N,1);
for i=1:N;
    nb_measures(i)=sum(ids==jtemp(i));
end;
ntotal=sum(nb_measures);
nb_design=max(nb_measures);
covariables=zeros(N,nb_covariables);
indicobs=zeros(N,nb_design);


jid=1;
jvarex1=jid+1;
jvarex2=jvarex1+nb_varex-1;
jy=jvarex2+1;
jcov1=jy+1;
jcov2=jcov1+nb_covariables-1;
i1=1;
indicobs=[];
y=S(:,jy);
id=zeros(size(S,1),1);
design=S(:,jvarex1:jvarex2);
for i=1:N;
    nbi=nb_measures(i);
    indicobs(i,1:nbi)=ones(1,nbi);
    covariables(i,:)=S(i1,jcov1:jcov2);
    id(i1:i1+nbi-1)=repmat(i,nbi,1);
    i1=i1+nbi;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grafsaem(PARA,kiter,hfig,Mstrm);

[niter, nb_param]=size(PARA);
n1=round(sqrt(nb_param));
n2=ceil(nb_param/n1);
for j=1:nb_param
    figure(hfig);
    subplot(n1,n2,j)
    %plot(PARA(1:kiter,j))
    semilogx(PARA(1:kiter,j))
    v=axis; v(1)=1; v(2)=niter; axis(v);
    set(gca,'fontsize',8)
    % title(Mstrm{j},'fontsize',10,'fontweight','bold')
    title(Mstrm{j},'fontsize',10)
    drawnow
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = normcdf(x,mu,sigma)
p = 0.5 * erfc(-x/sqrt(2));


