function [righths] = RHS(h)
Patm= 101325;a= 0.286561;%a=0.3432;
b= 0.046089;%b=0.0352;
p1= 0.45;p2= 2500;l=10;
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
Pvac=22700;sigma=Patm-Pvac;Vfoutlet=a*sigma^b;
houtlet=(p1*l)/(p2*Vfoutlet); %low value for map

%PDE Variables
N=20; %Number of Cheb points
%x = linspace(0,0.5,N);%(cos(pi*(0:N-1)/(N-1))'+1)./4;
%x=x(end:-1:1);
%h=h+[sin(0:2*pi/(N):2*pi)]/10000;
[x,D]=lagdif(N,2,320);
D2=D(:,:,2);
righths=Wfunc(p1,p2,C,n,l,a,b,mu,Patm,h).*(D2*h);