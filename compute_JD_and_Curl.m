function [f0, g0] = compute_JD_and_Curl(Phi1,Phi2,h)
[Phi1y, Phi1x] = gradient(Phi1,h);
[Phi2y, Phi2x] = gradient(Phi2,h);
f0 = Phi1x.*Phi2y-Phi1y.*Phi2x;
g0 = Phi2x-Phi1y;
end