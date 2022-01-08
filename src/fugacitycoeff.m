function [fugsmallcoeffi1,fugsmallcoeffi2] = fugacitycoeff(y1,...
  p,t,decide)
%y1 = 0;
y2 = 1-y1;
%ps = input("input pressure in kPa\n","S");
%if ps == 'd'
 % p = 95855;
%else
 % p = str2num(ps)*1000;
%end
%size(p)
%ts = input("input temperature\n","S");
%if ts == 'd'
 % t = 350;
%else
%  t = str2num(ts)
%end
%size(t)
p = p*1000;
R = 8.314;
pc1 = 51750;
pc2 = 45200;
tc1 = 536.8;
tc2 = 632.4;
pr1 = p/pc1;
pr2 = p/pc2;
tr1 = t/tc1;
tr2 = t/tc2;
acen1 = 0.622;
acen2 = 0.249;
%printf(["enter 1 to use Vander waals CEOS\n"...
%"enter 2 to use Redlich Kwong CEOS\nenter 3 to use Soave "...
%"Redlich Kwong CEOS\nenter 4 to use Peng Robinson CEOS\n"])
% ceos parameters
%decide = input("");
switch decide
  case 1  
  alpha1 = 1;
  alpha2 = 1;
  sigma = 0;
  epsilon = 0;
  omega = 1/8;
  psipar = 27/64;
  case 2  
  alpha1 = tr1^-0.5;
  alpha2 = tr2^-0.5;
  sigma = 1;
  epsilon = 0;
  omega = 0.08664;
  psipar = 0.42748;
  case 3
  alpha1...
  =( 1+ (0.480 + 1.574*acen1 - 0.176*acen1^2)*(1 -...
  tr1^(0.5)))^2;
  alpha2 =...
  (1+(0.480 + 1.574*acen2 - 0.176*acen2^2)*(1-tr2^(1/2)))^2;
  sigma = 1;
  epsilon = 0;
  omega = 0.08664;
  psipar = 0.42748;
  case 4
  alpha1...
  =( 1+ (0.37464 + 1.54226*acen1 - 0.26992*acen1^2 )*(1 -...
  tr1^(0.5)))^2;
  alpha2 =...
  (1+(0.37464+1.54226*acen2-0.26992*acen2^2)*(1-tr2^(1/2)))^2;
  sigma = 1+sqrt(2);
  epsilon = 1-sqrt(2);
  omega = 0.07779607;
  psipar = 0.45723553;
otherwise
error("please enter an integer from 1 to 4\n");
end
a1 = alpha1*psipar*R^2*tc1^2/pc1;
a2 = alpha2*psipar*R^2*tc2^2/pc2;
b1 = omega*R*tc1/pc1;
b2 = omega*R*tc2/pc2;
am = y1^2*a1+2*y1*y2*sqrt(a1*a2)+y2^2*a2;
bm = y1*b1+y2*b2;
a1bar = 2*y1*a1+ 2*y2*sqrt(a1*a2) - am;
a2bar = 2*y2*a2+ 2*y1*sqrt(a1*a2) - am;
b1bar = b1;
b2bar = b2;
betapar = bm*p/(R*t);
q = am/(bm*R*t);
q1bar = q*((2*y1*a1+ 2*y2*sqrt(a1*a2))/am - b1/bm);
q2bar = q*((2*y2*a2+ 2*y1*sqrt(a1*a2))/am - b2/bm);
A = (sigma+epsilon)*betapar-(betapar+1);
B =...
sigma*epsilon*betapar^2-(sigma+epsilon)*betapar*(betapar+1)...
+q*betapar;
C = -betapar^2*(q+sigma*epsilon*(betapar+1));
P = B - A^2/3;
Q = 2*A^3/27 - A*B/3 + C;
D = Q^2/4 + P^3/27;
if D > 0
  Z = nthroot((-Q/2+sqrt(D)),3) +nthroot((-Q/2-sqrt(D)),3)-A/3;
elseif D == 0
  Z = [-2*nthroot((Q/2),3)-A/3 nthroot((Q/2),3)-A/3];
elseif D<0
  ta = -Q/2*nthroot((-27/(P^3)),2);
  theta = acos(ta);
  r = nthroot((-(P^3)/27),2);
  Zs = [2*nthroot(r,3)*cos((theta)/3)-A/3,
  2*nthroot(r,3)*cos((2*pi+theta)/3)-A/3,
  2*nthroot(r,3)*cos((4*pi+theta)/3)-A/3];
  Z = [max(Zs) min(Zs)];
end
Z;
if decide == 1
  I = betapar/(Z+epsilon*betapar);
  else
  I =...
  1/(epsilon-sigma)*log((Z+epsilon*betapar)./(Z+sigma*betapar));
end

fugsmallcoeffi1 = exp(b1bar/bm*(Z-1)-log(Z-betapar)-q1bar*I);
fugsmallcoeffi2 = exp(b2bar/bm*(Z-1)-log(Z-betapar)-q2bar*I);
%fugacity1 = fugsmallcoeffi1*p*y1;
%fugacity2 = fugsmallcoeffi2*p*y2;
return;
endfunction