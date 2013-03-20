function K = Kfunc(C,n,p1,p2,l,h)
vf = vffunc(p1,p2,l,h);
K = C.*((1-vf).^(n+1))./(vf).^n;