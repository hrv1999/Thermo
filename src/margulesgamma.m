function [gammamarg] = margulesgamma(x1,y1,P,x,t)
pkg load optim;
%x1 = [0.055	0.129	0.212	0.313	0.430	0.520	0.638	0.749	0.872...
%];
%y1 = [0.350	0.511	0.599	0.650	0.697	0.726	0.759	0.813	0.883...
%];
%P = [46.66	59.46	69.06	76.53	81.19	84.32	88.66	90.86	92.19
%];
%t = 368.15;
x2 = 1 - x1;
y2 = 1 - y1;
an1 = [16.0692, 3448.66, -69.06];
an2 = [13.9926, 3295.12, -55.6];
p1sat = exp(an1(1)-an1(2)/(t+an1(3)));
p2sat = exp(an2(1)-an2(2)/(t+an2(3)));
gamma1exp = y1.*P./(x1*p1sat);
gamma2exp = y2.*P./(x2*p2sat);
gertexp = x1.*log(gamma1exp)+x2.*log(gamma2exp);
modelfunc = @(W,x) x.*(1-x).*(W(2)*x+W(1)*(1-x));
[W] = nonlin_curvefit(modelfunc,[1; 1],x1,gertexp);
gammamarg(1) = exp((1-x)^2*(W(1)+2*(W(2)-W(1))*x));
gammamarg(2) = exp(x^2*(W(2)+2*(W(1)-W(2))*(1-x)));
return;
endfunction