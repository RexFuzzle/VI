%Solving VI by PDEs
close all
clear all
clc
Patm= 101325;a= 0.286561;b= 0.046089;p1= 0.45;p2= 2500;l= 10;
C= 9.5e-11;n= 2.6;mu= 0.227;Pvac= 22700; 
% out   = map(in,inlow,inhigh,outlow,outhigh)
% outP  = Pfunc(Patm,a,b,p1,p2,l,h)
% vf    = vffunc(p1,p1,l,h)
% K     = Kfunc(C,n,p1,p2,l,h)
% dpdh  = dpdhfunc(p1,p2,l,a,b,h)
% dkdh  = dkdhfunc(p1,p2,l,C,n,h)
% W     = Wfunc(p1,p2,C,n,l,a,b,mu,Patm,h)
% [D,x] = cheb(N)
%BC 
Pinlet=90000;sigma=Patm-Pinlet;Vfinlet=a*sigma^b;
hinlet=(p1*l)/(p2*Vfinlet); %high value for map
Pvac=22500;sigma=Patm-Pvac;Vfoutlet=a*sigma^b;
houtlet=(p1*l)/(p2*Vfoutlet); %low value for map

%PDE Variables
N=20; %Number of Cheb points
h(1,:,1)=linspace(houtlet+10e-6,houtlet,N+1);
%h(1,1,1)=hinlet;
h(1,:,2)=h(1,:,1)+10e-8;
[D,x] = cheb(N);

D2=D^2;

titer=1;
niter=0;
err=1;
compiter(2)=1;
dx=0.01;

u=Kfunc(C,n,p1,p2,l,houtlet)/mu*(Pinlet-Pvac)/dx;

dt=dx/u

return

for t=dt:dt*2:dt*10
  titer++;
  niter=0;
  err(1)=1;
  err(2)=10;
  while abs(err(niter+2)-err(niter+1))>10e-8
    niter++;
    if niter>0
    h(niter+1,:,titer-1)= h(niter,:,titer-1);
    end
    temp=Wfunc(p1,p2,C,n,l,a,b,mu,Patm,h(niter,:,titer))';
    newH=(D2(2:end-1,2:end-1))\(((h(niter,2:end-1,titer)-h(compiter(titer),2:end-1,titer-1))/dt)./temp(2:end-1)')';
    h(compiter(titer),2:end-1,titer-1);
    newh=newH;%map(newH,1,-1,houtlet,hinlet);
    err(niter+2)=norm(h(niter,2:end-1,titer)-newh');
    h(niter+1,:,titer)=[0;newh;0]+((-x(end:-1:1)+1)/2*hinlet+(x(end:-1:1)+1)/2*houtlet);
  end
  %dx=
  %dpdx=(Pfunc(Patm,a,b,p1,p2,l,h())-Pfunc(Patm,a,b,p1,p2,l,h()))/dx;
  %dhdx=(h-h)/dx;
  %dhdt=(h-h)/dt;
  %u=-Kfunc(C,n,p1,p2,l,h())/mu*dPdx+((x()-xc())/h())*(Kfunc(C,n,p1,p2,l,h())/mu*dhdx*dpdx-dhdt);
  %dist=u*t
  %xstd=linspace(0,dist,N)
  %hstdx=interp1(map(x(end:-1:1),-1,1,0,dist),h(niter+1,:,titer),xstd)
  h(1,:,titer+1)=h(niter+1,:,titer);
  compiter(titer+1)=niter;
  titer
  niter
fflush(stdout);
end
hfinal(1,:)=h(1,:,1);
for temp=2:titer
hfinal(temp,:)=h(compiter(temp+1),:,temp);
end
figure(1)
plot(x(end:-1:1),hfinal');