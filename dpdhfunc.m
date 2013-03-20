function dpdh = dpdhfunc(p1,p2,l,a,b,h)
dpdh = (((p1*l)./(p2*h*a)).^(1/b))./(b*h);