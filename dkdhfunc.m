function dkdh = dkdhfunc(p1,p2,l,C,n,h)
dkdh=C.*(-(-p2.*h+p1*l)./p2./h).^n.*(p1*l./p2./h).^(-n).*(p1*l+n*p2.*h)./p2./h.^2;
