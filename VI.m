%Solving VI by PDEs
close all
clear all
clc
Patm= 101325;a= 0.286561;%a=0.3432;
b= 0.046089;%b=0.0352;
p1= 0.45;p2= 2500;l=10;
C= 9.5e-11;n= 2.6;mu= 0.227;Pvac= 22700; 
length=0.5;
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
N=100;
tsteps=10e9;
h=zeros(100,100);
h(2,1)=houtlet;
h(1,1)=houtlet;
h(1,2)=hinlet;
h(2,2)=houtlet;
dx=length/N;
%h(1:10,10)=linspace(hinlet,houtlet,10)';
dtp=0;
for ff=3:100
x=linspace(0,dx*ff,ff);
D=poldif(x,2);
D2=D(:,:,2);
%D2=D2(2:end-1,2:end-1);
h(:,ff)=h(:,ff-1);
h(ff,ff)=houtlet;
u=Kfunc(C,n,p1,p2,l,houtlet)/mu*((Pinlet-Pvac)/(dx*ff))
dt=((dx*ff)/u)-dtp
ht=h(1:ff,ff);
titer=1;
for dtit=dt/tsteps:dt/tsteps*1000000:dt
  titer++;
  temp=Wfunc(p1,p2,C,n,l,a,b,mu,Patm,ht(:,titer-1))';
  ht(1:end,titer)=(dt/tsteps*temp(1:end)'.*(D2*ht(1:end,titer-1)))+ht(1:end,titer-1);
  ht(1,titer)=hinlet;
  ht(end,titer)=houtlet;
  fflush(stdout);
%plot(ht(:,end))
%drawnow
titer;
end
h(1:ff,ff)=ht(:,titer);
%h(2,3)=0.0039;
h;
dtp=dt;
end

return
figure(1)
%plot(x(end:-1:1),h(1:end,:)');
hplot=[1:100:length(h)];
mesh(x(end:-1:1),1:100:length(h),h(hplot,:))
view(45,45);