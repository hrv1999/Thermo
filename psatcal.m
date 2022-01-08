function [psat] = psatcal(t)
an1 = [16.0692, 3448.66, -69.06];
an2 = [13.9926, 3295.12, -55.6];
psat(1) = exp(an1(1)-an1(2)/(t+an1(3)));
psat(2) = exp(an2(1)-an2(2)/(t+an2(3)));
  return;
  endfunction