function [ tau ] = transition_tilt( t,psi1,s,psi2 )
%TRANSITION_TILT Summary of this function goes here
%   Detailed explanation goes here

cos_2 = (cos(psi1-psi2))^2;
g = ( (t/s)^2 + 1 )*cos_2 + (1/s^2 + t^2)*(1-cos_2);
G = (s/t)*g/2;
tau = G + sqrt(G^2 -1);

end

