% function that returs capital phi k
function [capphik] = capitalphik(y1,p,t)
  %an represents antoine constants
  an1 = [16.0692, 3448.66, -69.06];
  an2 = [13.9926, 3295.12, -55.6];
  p1sat = exp(an1(1)-an1(2)/(t+an1(3)));
  p2sat = exp(an2(1)-an2(2)/(t+an2(3)));
  capphik(1)=fugacitycoeff(y1,p,t,4)./fugacitycoeff(1,p1sat,t,4);
  capphik(2)=fugacitycoeff(1-y1,p,t,4)./fugacitycoeff(0,p2sat,t,4);
  return
  endfunction