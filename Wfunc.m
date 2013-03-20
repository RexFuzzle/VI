function W = Wfunc(p1,p2,C,n,l,a,b,mu,Patm,h)
W = (1/mu)*(dkdhfunc(p1,p2,l,C,n,h)+(Kfunc(C,n,p1,p2,l,h).*(1+h))).*dpdhfunc(p1,p2,l,a,b,h);
