%Non linear finite difference example
clear all;close all;clc;
n=16;eps=1;dt=0.0001;
[D,x]=cheb(n);D2=D^2;
%D2(end,:)=D(end,:)-1;
%D2(1,:)=0;
%D2(:,1)=0;
u(:,1)=linspace(0,2,n+1);
u(end,1)=1*(x(end-1)-x(end))+u(end-1,1);
tsiter=0;
for t=dt:dt:0.0009
  tsiter++;
  u(1:n+1,tsiter+1)=(-eps*dt*exp(-u(1:n+1,tsiter)).*(D2*u(1:n+1,tsiter)))+u(1:n+1,tsiter);
  u(1,tsiter+1)=0;
  u(end,tsiter+1)=2;
  %u(end,tsiter+1)=1*(x(end-1)-x(end))+u(end-1,tsiter);
end
surf(1:length(u),x(end:-1:1),u(:,1:end))
%plot(linspace(0,1,n),u_new)