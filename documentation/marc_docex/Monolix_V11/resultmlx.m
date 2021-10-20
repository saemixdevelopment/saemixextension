function  resultmlx(varargin);
%RESULTMLX  graphical display of the results, called by MONOLIX
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

if isstruct(varargin{1})==0
    if  strcmp(varargin{1},'graf7mlx')==1
        graf7mlx(varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
    elseif strcmp(varargin{1},'graf2mlx')==1
        graf2mlx(varargin{2});
    elseif strcmp(varargin{1},'graf4mlx')==1
        graf4mlx(varargin{2});
    elseif strcmp(varargin{1},'graf6mlx')==1
        graf6mlx(varargin{2},varargin{3},varargin{4});
    end;
else
    s=varargin{1};
    res=s.resultats;
    if nargin<3
        close all
        ngraf=7;
        hfig=zeros(ngraf,1);
        for k=1:ngraf
            hfig(k)=figure('visible','off');
        end;
        addmenumlx(hfig);
    else
        hfig=varargin{3};
    end;
    if nargin<2
        K=500;
    else
        K=varargin{2};
    end;

    load(res);
    Gamma2_phi1=cov_random(i1,i1);
    vprior_phi1=diag(Gamma2_phi1)';

    n1max=3;
    n2max=4;

    indio=find(io');
    [N,nt_max]=size(io);
    nphi1=length(i1);
    nphi0=length(i0);
    nphi=nphi1+nphi0;

    mprior_phi1=mprior_phi(:,i1,:);
    mprior_phi0=mprior_phi(:,i0,:);
    [fprior,gprior]=feval(fonction,mprior_phi,design,id);
    [fpost,gpost]=feval(fonction,mpost_phi,design,id);
    nmax=n1max*n2max;
    %Nmax=ceil(N/nmax);
    nx=size(design,2);
    temp=zeros(N*nt_max,1);
    temp(indio)=y;
    y=reshape(temp,nt_max,N)';
    temp(indio)=fprior;
    fprior=reshape(temp,nt_max,N)';
    temp(indio)=fpost;
    fpost=reshape(temp,nt_max,N)';
    t=zeros(N,nt_max,nx);
    for ix=1:nx
        temp(indio)=design(:,ix);
        t(:,:,ix)=reshape(temp,nt_max,N)';
    end;

    if nx>1
        nb_varex=size(t,3);
        str1=cell(nb_varex,1);
        str2=cell(nb_varex+1,1);
        str2{1}='none ';
        for k=1:nb_varex
            str1{k}=['var ' num2str(k)];
            str2{k+1}=['var ' num2str(k)];
        end;
        nvar=menu('regression variable:',str1);
        xvar=t(:,:,nvar);
        %ncat=menu('categorical variable:',str2)-1;
        ncat=0;
        if ncat>0
            xcat=t(:,:,ncat);
            ucat=unique(xcat);
            disp(ucat(:)');
            k=menu('value of the categorical variable:',num2cell(ucat));
            vcat=ucat(k);
        else
            xcat=ones(size(y(:,:,1)));
            vcat=1;
        end;
    else
        xvar=t;
        xcat=ones(size(y(:,:,1)));
        vcat=1;
        nvar=1;
        ncat=0;
    end;

    n1=min(fix(sqrt(N)),n1max);
    n2=min(ceil(N/n1),n2max);

    vaxisy=zeros(4,1);
    v1=min(min(xvar(1:n1*n2,:)));
    v2=max(max(xvar(1:n1*n2,:)));
    dx=v2-v1;
    vaxisy(1)=v1-dx/10;
    vaxisy(2)=v2+dx/10;
    i12=1:n1*n2;
    y12=y(i12,:);
    jy12=find( (io(i12,:)==1) & (xcat(i12,:)==vcat) );
    v3=min(min(y12(jy12)));
    v4=max(max(y12(jy12)));
    dy=v4-v3;
    vaxisy(3)=v3-dy/10;
    vaxisy(4)=v4+dy/10;
    jo=find ((io==1) & (xcat==vcat));
    %
    nphi1=length(i1);
    nphi0=length(i0);
    nphi=nphi1+nphi0;
    mtild_phi1=repmat(mprior_phi1,[1,1,K]);
    vtild_phi1=repmat(vprior_phi1,[N,1,K]);
    phi=zeros(N,nphi,K);
    phi(:,i1,:)=mtild_phi1+sqrt(vtild_phi1).*repmat(randn(1,nphi1,K),[N,1,1]);
    phi(:,i0,:)=repmat(mpost_phi(:,i0),[1,1,K]);
    [fsim,gsim]=feval(fonction,phi,repmat(design,[1,1,K]),id);
    ysim=fsim+(gsim.*randn(size(fsim)))*sqrt(sigma2);
    ytemp=zeros(N,nt_max,K);
    for ix=1:K
        temp(indio)=ysim(:,ix);
        ytemp(:,:,ix)=reshape(temp,nt_max,N)';
    end;
    ysim=ytemp;
    ymean=mean(ysim,3);
    sey=std(ysim,0,3);
    %
    figure(hfig(1))
    for i=1:N
        ji=find( (io(i,:)==1) & (xcat(i,:)==vcat) );
        ti=xvar(i,ji);
        yi=y(i,ji);
        plot(ti,yi,'o-')
        hold on
    end;
    vax=axis; vax(1:2)=vaxisy(1:2);
    axis(vax);
    xlabel('x');
    ylabel('observations');
    hold off;
    %
    figure(hfig(2))
    save temp2 n1 n2  y io fprior fpost ymean sey xvar xcat vcat
    i0=1;
    graf2mlx(i0)

    %
    figure(hfig(3))
    v1=min([y(jo) ; fprior(jo) ; fpost(jo)]);
    v2=max([y(jo) ; fprior(jo) ; fpost(jo)]);
    v=[v1 v2 v1 v2];
    subplot(221)
    plot(y(jo),fprior(jo),'.')
   % loglog(y(jo),fprior(jo),'.')
    axis(v);
    line([v1 v2],[v1 v2],'color','red')
    xlabel('observations')
    ylabel('predictions')
    title('Using the population parameters')
    subplot(222)
    plot(y(jo),fpost(jo),'.')
    %loglog(y(jo),fpost(jo),'.')
    axis(v);
    line([v3 v2],[v3 v2],'color','red')
    xlabel('observations')
    ylabel('predictions')
    title('Using the individual parameters')
    subplot(223)
    plot(xvar(jo),y(jo)-fprior(jo),'.')
    v=axis;
    line([v(1) v(2)],[0 0],'color','red')
    grid
    xlabel('x')
    ylabel('residuals')
    subplot(224)
    plot(xvar(jo),y(jo)-fpost(jo),'.')
    v=axis;
    line([v(1) v(2)],[0 0],'color','red')
    grid
    ylabel('residuals')
    xlabel('x')
    %
    figure(hfig(4))
    i0=1;
    graf4mlx(i0)
    %
    figure(hfig(5))
    plot(xvar(jo),(y(jo)-ymean(jo))./sey(jo),'.')
    v=axis;
    line([v(1) v(2)],[0 0],'color','red')
    grid
    ylabel('residuals')
    xlabel('x')
    %
    ncovariables=size(covariables,2);
    if ncovariables >0
        figure(hfig(6))
        clf
        addmenumlx(gcf);
        i2=min(i1);
               i1=1;
 P1=[.01 .8 .08 .04];
        P2=[.01 .6 .08 .2];
        P3=[.01 .4 .08 .04];
        P4=[.01 .2 .08 .2];
        for k=1:ncovariables
            stlist1{k}=['c',num2str(k)];
        end;
        for k=1:nphi
            stlist2{k}=['phi',num2str(k)];
        end;
         uicontrol('unit','normalized','style','text','string','Covariate','position',P1,'fontsize',8,'FontUnits','normalized');
        h1=uicontrol('unit','normalized','style','listbox','string',stlist1,'position',P2,'fontsize',8,'FontUnits','normalized','value',i1);
        uicontrol('unit','normalized','style','text','string','Parameter','position',P3,'fontsize',8,'FontUnits','normalized');
        h2=uicontrol('unit','normalized','style','listbox','string',stlist2,'position',P4,'fontsize',8,'FontUnits','normalized','value',i2);
%        uicontrol('unit','normalized','style','text','string','Covariate','position',P1,'fontweight','bold','fontsize',8);
%         h1=uicontrol('unit','normalized','style','listbox','string',stlist1,'position',P2,'fontweight','bold','fontsize',8,'value',i1);
%         uicontrol('unit','normalized','style','text','string','Parameter','position',P3,'fontweight','bold','fontsize',8);
%         h2=uicontrol('unit','normalized','style','listbox','string',stlist2,'position',P4,'fontweight','bold','fontsize',8,'value',i2);
        %strclk1=['i1=str2num(get(gco,''value''); h2=get(gco,''userdata''); i2=str2num(get(h2,''value''); resultmlx(''graf6mlx'',res,i1,i2);'];
        strclk=['h=get(gco,''userdata''); i2=get(h{2},''value'');  i1=get(h{1},''value''); resultmlx(''graf6mlx'',h{3},i1,i2);'];
        h{1}=h1;
        h{2}=h2;
        h{3}=res;
        set(h1,'callback', strclk,'userdata',h);
        set(h2,'callback', strclk,'userdata',h);
        graf6mlx(res,i1,i2)
        set(hfig(6),'name','  Covariates');
    end;
    %
    graf7mlx(res,nvar,ncat,vcat,0,hfig(7));

    set(hfig(1),'name','  Spaghetti plot');
    set(hfig(2),'name','  Indivudal fits');
    set(hfig(3),'name','  Predictions vs observations');
    set(hfig(4),'name','  Population distributions');
    set(hfig(5),'name','  Population residuals');
    set(hfig(7),'name','  Data and model');

    vermat=version;
    if str2num(vermat(1))>=7
        set(hfig,'WindowStyle','docked')
    end;

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function graf2mlx(i0)
clf
addmenumlx(gcf);
P1=[.01 .6 .07 .04];
uicontrol('unit','normalized','style','text','string','first ID','position',P1,'fontsize',8,'FontUnits','normalized');
P2=[.01 .56 .07 .04];
h=uicontrol('unit','normalized','style','edit','string',num2str(i0),'position',P2,'fontweight','bold','fontsize',8,'FontUnits','normalized');
set(h,'callback','i0=str2num(get(gco,''string'')); resultmlx(''graf2mlx'',i0);');
load temp2;
N=size(y,1);
ndisp=min(n1*n2,N-i0+1);
for ip=1:ndisp
    subplot(n1,n2,ip);
    i=ip+i0-1;
    ji=find( (io(i,:)==1) & (xcat(i,:)==vcat) );
    ti=xvar(i,ji);
    yi=y(i,ji);
    f1i=fprior(i,ji);
    f2i=fpost(i,ji);
    if length(ji)<=1
        plot(ti,yi,'+',ti,f1i,'ro-',ti,f2i,'go-')
    else
        plot(ti,yi,'+',ti,f1i,'r-',ti,f2i,'g-')
    end;
    d=max(ti)-min(ti);
    if d>0
        vax=axis; vax(1:2)=[min(ti)-d/10; max(ti)+d/10];
        axis(vax);
    end;
    title(i)
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function graf4mlx(i0)
clf
addmenumlx(gcf);
P1=[.01 .6 .07 .04];
uicontrol('unit','normalized','style','text','string','first ID','position',P1,'fontsize',8,'FontUnits','normalized');
P2=[.01 .56 .07 .04];
h=uicontrol('unit','normalized','style','edit','string',num2str(i0),'position',P2,'fontweight','bold','fontsize',8,'FontUnits','normalized');
set(h,'callback','i0=str2num(get(gco,''string'')); resultmlx(''graf4mlx'',i0);');
load temp2
N=size(y,1);
ndisp=min(n1*n2,N-i0+1);
for ip=1:ndisp
    subplot(n1,n2,ip);
    i=ip+i0-1;
    ji=find( (io(i,:)==1) & (xcat(i,:)==vcat) );
    ti=xvar(i,ji);
    yi=y(i,ji);
    f0i=ymean(i,ji);
    f1i=ymean(i,ji)-sey(i,ji);
    f2i=ymean(i,ji)+sey(i,ji);
    if length(ji)<=1
        plot(ti,yi,'+',ti,f0i,'ro',ti,[f1i' f2i'],'r+')
    else
        plot(ti,yi,'+',ti,f0i,'r-',ti,[f1i' f2i'],'r:')
    end;
    d=max(ti)-min(ti);
    if d>0
        vax=axis; vax(1:2)=[min(ti)-d/10; max(ti)+d/10];
        axis(vax);
    end;
    title(i)
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function graf6mlx(res,ik1,ik2)
load(res)
N=size(io,1);
nphi=size(mpost_phi,2);
x=covariables(:,ik1);
y1=mpost_phi(:,ik2);
y2=mpost_phi(:,ik2)-mprior_phi(:,ik2);
C=cov([x,y1,y2]);
if C(1)*C(5)==0
    r1=NaN;
else
r1=C(2)/sqrt(C(1)*C(5));
end;
if C(1)*C(9)==0
    r2=NaN;
else
r2=C(3)/sqrt(C(1)*C(9));
end;
a1=C(2)/C(1);
b1=mean(y1)-a1*mean(x);
a2=C(3)/C(1);
b2=mean(y2)-a2*mean(x);
subplot(211)
plot(x,y1,'.')
v=axis;
plot(x,y1,'.',[v(1) v(2)],[a1*v(1)+b1 a1*v(2)+b1],'-')
axis(v);
title(['correlation coefficient: ', num2str(r1,'%5.4f')])
xlabel(['c',num2str(ik1)])
ylabel(['E(phi',num2str(ik2),'|y)'])
set(gca,'position',[0.2 0.61 0.7 0.33])
subplot(212)
plot(x,y2,'.')
v=axis;
plot(x,y2,'.',[v(1) v(2)],[a2*v(1)+b2 a2*v(2)+b2],'-')
axis(v);
title(['correlation coefficient: ', num2str(r2,'%5.4f')])
xlabel(['c',num2str(ik1)])
ylabel(['E(b',num2str(ik2),'|y)'])
set(gca,'position',[0.2 0.09 0.7 0.33])
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function graf7mlx(res,nvar,ncat,vcat,iop,h);

if exist('iop')==0 | iop==0
    nx=200;
    load(res)
    nb_varex=size(design,2);

    if nb_varex>1
        if exist('nvar')==0
            str1=cell(nb_varex,1);
            str2=cell(nb_varex+1,1);
            str2{1}='none ';
            for k=1:nb_varex
                str1{k}=['var ' num2str(k)];
                str2{k+1}=['var ' num2str(k)];
            end;
            nvar=menu('regression variable:',str1);
            ncat=menu('categorical variable:',str2)-1;
            if ncat>0
                xcat=design(:,ncat);
                ucat=unique(xcat);
                disp(ucat(:)');
                k=menu('value of the categorical variable:',num2cell(ucat));
                vcat=ucat(k);
            else
                xcat=ones(size(design(:,1)));
                vcat=1;
            end;
        else
            xcat=ones(size(design(:,1)));
            vcat=1;
        end;
    else
        nvar=1;
        xcat=ones(size(design(:,1)));
        vcat=1;
    end;
    j=find(  xcat==vcat );
    y=y(j,:);
    %mpost_phi=mpost_phi(j,:);
    mmprior_phi=mean(mprior_phi,1);
    design=design(j,:);
    xvar=design(:,nvar);
    T=cumsum([0;nb_measures(1:N-1)])+1;
    if exist('h')==0
        h=figure;
        addmenumlx(h);
    else
        figure(h);
    end;
    plot(xvar,y);
    v=axis;
    x=linspace(v(1),v(2),nx)';
    tm=repmat(mean(design,1),nx,1);
    tm(:,nvar)=x;
    [f,g]=feval(fonction,repmat(mmprior_phi,nx,1),tm,1:nx);
    pl1=plot(xvar,y,'.',x,f);
    xlabel('x');
    ylabel('y');
    title('click on a point to see the individual fit')
    legend('observed data','population model')
    save temp7  f tm design mmprior_phi mpost_phi x y id fonction xvar nx nvar T
    set(pl1(1),'ButtonDownFcn','resultmlx(''graf7mlx'',[],[],[],[],1)','UserData',gca);
    set(gca,'ButtonDownFcn','resultmlx(''graf7mlx'',[],[],[],[],2)','UserData',pl1(1));
    %  set(pl1(1),'ButtonDownFcn','graf7mlx([],[],[],[],1);','UserData',gca);
    %   set(gca,'ButtonDownFcn','graf7mlx([],[],[],[],2);','UserData',pl1(1));
    drawnow
    %%%%%%%%%%%%%%%%%%%%%%
else
    load temp7
    if iop==1
        h=gcbo;
        a=get(h,'UserData');
    else
        a=gcbo;
        h=get(a,'UserData');
    end;
    xdata=get(h,'XData');
    ydata=get(h,'YData');
    PositionClick=get(a,'CurrentPoint');
    X=PositionClick(1,1);
    Y=PositionClick(1,2);
    v=axis;
    dx=(xdata-X)/(v(2)-v(1));
    dy=(ydata-Y)/(v(4)-v(3));
    distance=sqrt(dx.^2+dy.^2);
    [DistanceMin,index]=min(distance);
    j=id(index);
    tj=repmat(design(T(j),:),nx,1);
    tj(:,nvar)=x;
    [fj,g]=feval(fonction,repmat(mpost_phi(j,:),nx,1),tj,1:nx);
    %[fj,g]=saquinavir_funct(mpost_phi(j,:),x,1);
    idj=find(id==j);
    pl1=plot(xvar,y,'.',x,f,'-',x,fj,'r-',xvar(idj),y(idj),'or');
    title('click on a point to see the individual fit')
    xlabel('x');
    ylabel('y');
    legend('observed data','population model','individual model')
    set(pl1(1),'ButtonDownFcn','resultmlx(''graf7mlx'',[],[],[],[],1)','UserData',gca);
    set(gca,'ButtonDownFcn','resultmlx(''graf7mlx'',[],[],[],[],2)','UserData',pl1(1));
    %   set(pl1(1),'ButtonDownFcn','graf7mlx([],[],[],[],1);','UserData',gca);
    %  set(gca,'ButtonDownFcn','graf7mlx([],[],[],[],2);','UserData',pl1(1));
end;


