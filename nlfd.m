%Non linear finite difference example
clear all;close all;clc;
n=6;h=1/(n+1);eps=1;
row=[2 -1 zeros(1,n-2)];
A=toeplitz(row,row);
A(1,:)=0;
A(end,end)=h;A(end,end-1)=-h;
A=A.*(h^-2).*eps;A(1,1)=1;
u_old=(linspace(0,0.1,n))';
diff=1;
iter=0;
while diff>10e-6
iter++;
b=exp(u_old)';b(1)=0;
b(end)=1+h/2*b(end);
u_new=A\b';
diff=norm(u_new-u_old);
u_old=u_new;
end
plot(linspace(0,1,n),u_new)