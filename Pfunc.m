function outP = Pfunc(Patm,a,b,p1,p2,l,h)
vf = vffunc(p1,p2,l,h);
outP = Patm - (vf./a).^(1/b);