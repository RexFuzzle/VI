%Non linear finite difference example
clear all;close all;clc;
n=30;h=1/(n+1);eps=1;dt=0.001;
row=[2 -1 zeros(1,n-2)];
A=toeplitz(row,row);A(1,:)=0;
A(end,end)=0;A(end,end-1)=0;
A=A.*(h^-2).*eps;A(1,1)=1;
u(:,1)=ones(n,1)+1;
tsiter=0;
for t=dt:dt:3
  tsiter++;
  u(:,tsiter+1)=(-eps*dt*exp(-u(:,tsiter)).*(A*u(:,tsiter)))+u(:,tsiter);
  u(1,tsiter+1)=0;
  u(end,tsiter+1)=h+u(end-1,tsiter+1);
end
surf(u);
%plot(linspace(0,1,n),u_new)