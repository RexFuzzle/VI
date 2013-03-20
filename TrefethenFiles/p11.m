% p11.m - Chebyshev differentation of a smooth function

  xx = 0:.001:1; uu = exp(xx).*sin(5*xx); clf
  for N = [10 20]
    [D,x] = cheb(N); u = exp(x).*sin(5*x);
     subplot(2,2,N/10)
%subplot('position',[.15 .66-.4*(N==20) .31 .28])
      plot(x,u,'.*','markersize',2), grid on
      line(xx,uu)
      title(['u(x),  N=' int2str(N)])
    error = D*u - exp(x).*(sin(5*x)+5*cos(5*x));
      subplot(2,2,N/10+2)
%subplot('position',[.55 .66-.4*(N==20) .31 .28])
      plot(x,error,'*','markersize',2), grid on
      line(x,error)
      title(['  error in u''(x),  N=' int2str(N)])
  end
