function Coeffs = PARdGroebnerWsCoefs(k,WT,E2T)
% this function calculates the coefficient of the groebner basis for Ws

a1 = k(1); d1 = k(2); k1 = k(3);
a2 = k(4); d2 = k(5); k2 = k(6);
a3 = k(7); d3 = k(8); k3 = k(9);
K1 = (d1+k1)/a1; K2 = (d2+k2)/a2; K3 = (d3+k3)/a3;
Wt = WT; E2 = E2T;

Coeffs = zeros(1,10);

% Ws^5
Coeffs(1) = (-k1*K3*k3 + 2*k1*K1*k3 - k3^2*K1);
% Ws^4 E1T
Coeffs(2) = K3*(2*k1^2 - k1*k3);
% Ws^4
Coeffs(3) = (-2*k1*k2*E2*K3 + 2*k1*k3*Wt*K3 - 2*k1*k3*E2*K3 - 2*k1*k3*K2*K3 - k2*k3*E2*K3 + 4*k1*k2*E2*K1 - 2*k1*k3*Wt*K1 + 2*k1*k3*E2*K1 + k1*k3*K1*K3 + 4*k1*k3*K1*K2 - 2*k2*k3*E2*K1 + k3^2*Wt*K1 - k3^2*E2*K1 - 2*k3^2*K1*K2);
% Ws^3 E1T
Coeffs(4) = K3*(-2*k1^2*Wt + k1*k3*Wt + 2*k1^2*E2 +2*k1*k2*E2 - k1*k3*E2 - k2*k3*E2 + k1^2*K3 + 4*k1^2*K2 - 2*k1*k3*K2);
% Ws^3 1e-16% error
Coeffs(5) = (2*k1*k2*Wt*E2*K3 - 2*k1*k2*E2^2*K3 - k1*k2*E2*K3^2 - k1*k3*Wt^2*K3 + 2*k1*k3*Wt*E2*K3 - k1*k3*E2^2*K3 - 2*k1*k2*E2*K2*K3 + 4*k1*k3*Wt*K2*K3 - 2*k1*k3*E2*K2*K3 -k1*k3*K2^2*K3 - 2*k2^2*E2^2*K3 + k2*k3*Wt*E2*K3 - k2*k3*E2^2*K3 - k2*k3*E2*K2*K3 + 4*k1*k2*E2*K1*K3 - k1*k3*Wt*K1*K3 + k1*k3*E2*K1*K3 + 4*k1*k2*E2*K1*K2 - 4*k1*k3*Wt*K1*K2 + 2*k1*k3*E2*K1*K2 + 2*k1*k3*K1*K2*K3 + 2*k1*k3*K1*K2^2 - k2*k3*E2*K1*K3 - 2*k2*k3*E2*K1*K2 + 2*k3^2*Wt*K1*K2 - k3^2*E2*K1*K2 - k3^2*K1*K2^2);
% Ws^2 E1T 1e-16% error
Coeffs(6) = K3*(-k1^2*Wt*K3 + k1^2*E2*K3 - 4*k1^2*Wt*K2 + 2*k1^2*E2*K2 + 2*k1^2*K2*K3 + 2*k1^2*K2^2 + k1*k2*E2*K3 + 2*k1*k2*E2*K2 + 2*k1*k3*Wt*K2 - k1*k3*E2*K2 - k1*k3*K2^2 - k2*k3*E2*K2);
% Ws^2 
Coeffs(7) = (k1*k2*Wt*E2*K3^2 - k1*k2*E2^2*K3^2 + 2*k1*k2*Wt*E2*K2*K3 - k1*k2*E2*K2*K3^2 - 2*k1*k3*Wt^2*K2*K3 + 2*k1*k3*Wt*E2*K2*K3 + 2*k1*k3*Wt*K2^2*K3 - k2^2*E2^2*K3^2 + k2*k3*Wt*E2*K2*K3 + k1*k2*E2*K1*K3^2 + 4*k1*k2*E2*K1*K2*K3 - 2*k1*k3*Wt*K1*K2*K3 + k1*k3*E2*K1*K2*K3 - 2*k1*k3*Wt*K1*K2^2 + k1*k3*K1*K2^2*K3 - k2*k3*E2*K1*K2*K3 + k3^2*Wt*K1*K2^2);
% Ws E1T 3e-16% error
Coeffs(8) = K3*(-2*k1^2*Wt*K2*K3 - 2*k1^2*Wt*K2^2 + k3*k1*Wt*K2^2 + k1^2*E2*K2*K3 + k2*k1*E2*K2*K3 + k1^2*K2^2*K3);
% Ws
Coeffs(9) = k1*K2*K3*(k2*Wt*E2*K3 - k3*Wt^2*K2 + k2*E2*K1*K3 - k3*Wt*K1*K2);
% E1T 1e-16% error
Coeffs(10) = -k1^2*K3^2*K2^2*Wt;

end