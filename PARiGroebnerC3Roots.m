function Coefficients = PARiGroebnerC3Roots(E2T,WT,k)
% this function calculates the coefficient of the groebner basis for C3

a1 = k(1); d1 = k(2); k1 = k(3);
a2 = k(4); d2 = k(5); k2 = k(6);
a3 = k(7); d3 = k(8);
a4 = k(9); d4 = k(10); k4 = k(11);

Coefficients = zeros(1,31);

% 1 C3 as a function of C4,Ws
Coefficients(1) = (-a1*d1*k1*a2*d3^2*a4*d4-a1*d1*k1*a2*d3^2*a4*k4+a1*d1*k1*d2*a3*d3*a4*d4+a1*d1*k1*d2*a3*d3*a4*k4+a1*d1*k1*k2*a3*d3*a4*d4+a1*d1*k1*k2*a3*d3*a4*k4-a1*k1^2*a2*d3^2*a4*d4-a1*k1^2*a2*d3^2*a4*k4+a1*k1^2*d2*a3*d3*a4*d4+a1*k1^2*d2*a3*d3*a4*k4+a1*k1^2*k2*a3*d3*a4*d4+a1*k1^2*k2*a3*d3*a4*k4+d1^2*a2*d3^2*a4^2*k4-d1^2*d2*a3*d3*a4^2*k4-d1^2*k2*a3*d3*a4^2*k4+2*d1*k1*a2*d3^2*a4^2*k4-2*d1*k1*d2*a3*d3*a4^2*k4-2*d1*k1*k2*a3*d3*a4^2*k4+k1^2*a2*d3^2*a4^2*k4-k1^2*d2*a3*d3*a4^2*k4-k1^2*k2*a3*d3*a4^2*k4);
% *C3*E1T+
Coefficients(2) = (a1*d1*a2*k2*d3^2*a4*d4*E2T+a1*d1*a2*k2*d3^2*a4*k4*E2T+a1*k1*a2*k2*d3^2*a4*d4*E2T+a1*k1*a2*k2*d3^2*a4*k4*E2T-d1^2*a2*k2*d3^2*a4^2*E2T-2*d1*k1*a2*k2*d3^2*a4^2*E2T-k1^2*a2*k2*d3^2*a4^2*E2T);
% *C3+
Coefficients(3) = (-d1^2*a2*a3^2*a4^2*k4-2*d1*k1*a2*a3^2*a4^2*k4-k1^2*a2*a3^2*a4^2*k4);
% *C4^2*Ws^2+
Coefficients(4) = (-a1*d1*k1*a2*a3*d3*a4*d4-a1*d1*k1*a2*a3*d3*a4*k4-a1*k1^2*a2*a3*d3*a4*d4-a1*k1^2*a2*a3*d3*a4*k4-d1^2*d2*a3^2*a4^2*k4-d1^2*k2*a3^2*a4^2*k4-2*d1*k1*d2*a3^2*a4^2*k4-2*d1*k1*k2*a3^2*a4^2*k4-k1^2*d2*a3^2*a4^2*k4-k1^2*k2*a3^2*a4^2*k4);
% *C4^2*Ws+
Coefficients(5) = (-a1*d1*k1*d2*a3*d3*a4*d4-a1*d1*k1*d2*a3*d3*a4*k4-a1*d1*k1*k2*a3*d3*a4*d4-a1*d1*k1*k2*a3*d3*a4*k4-a1*k1^2*d2*a3*d3*a4*d4-a1*k1^2*d2*a3*d3*a4*k4-a1*k1^2*k2*a3*d3*a4*d4-a1*k1^2*k2*a3*d3*a4*k4);
% *C4^2+
Coefficients(6) = (-d1^2*a2*a3^2*a4^2*k4-2*d1*k1*a2*a3^2*a4^2*k4-k1^2*a2*a3^2*a4^2*k4);
% *C4*Ws^3+
Coefficients(7) = (-a1*d1*k1*a2*a3*d3*a4*d4-a1*d1*k1*a2*a3*d3*a4*k4-a1*k1^2*a2*a3*d3*a4*d4-a1*k1^2*a2*a3*d3*a4*k4+d1^2*a2*a3^2*a4^2*k4*WT-d1^2*a2*a3^2*a4^2*k4*E2T-d1^2*d2*a3^2*a4^2*k4-d1^2*k2*a3^2*a4^2*k4+2*d1*k1*a2*a3^2*a4^2*k4*WT-2*d1*k1*a2*a3^2*a4^2*k4*E2T-2*d1*k1*d2*a3^2*a4^2*k4-2*d1*k1*k2*a3^2*a4^2*k4+k1^2*a2*a3^2*a4^2*k4*WT-k1^2*a2*a3^2*a4^2*k4*E2T-k1^2*d2*a3^2*a4^2*k4-k1^2*k2*a3^2*a4^2*k4);
% *C4*Ws^2+
Coefficients(8) = (-a1*d1*a2*a3*d3*a4*d4*k4-a1*d1*a2*a3*d3*a4*k4^2-a1*k1*a2*a3*d3*a4*d4*k4-a1*k1*a2*a3*d3*a4*k4^2+d1^2*a2*a3*d3*a4^2*k4+2*d1*k1*a2*a3*d3*a4^2*k4+k1^2*a2*a3*d3*a4^2*k4);
% *C4*Ws*E1T+
Coefficients(9) = (a1*d1*k1*a2*a3*d3*a4*d4*WT-a1*d1*k1*a2*a3*d3*a4*d4*E2T+a1*d1*k1*a2*a3*d3*a4*k4*WT-a1*d1*k1*a2*a3*d3*a4*k4*E2T-a1*d1*k1*d2*a3*d3*a4*d4-a1*d1*k1*d2*a3*d3*a4*k4-a1*d1*k1*k2*a3*d3*a4*d4-a1*d1*k1*k2*a3*d3*a4*k4-a1*d1*a2*k2*a3*d3*a4*d4*E2T-a1*d1*a2*k2*a3*d3*a4*k4*E2T-a1*d1*a2*a3*d3*a4*d4*k4*E2T-a1*d1*a2*a3*d3*a4*k4^2*E2T+a1*k1^2*a2*a3*d3*a4*d4*WT-a1*k1^2*a2*a3*d3*a4*d4*E2T+a1*k1^2*a2*a3*d3*a4*k4*WT-a1*k1^2*a2*a3*d3*a4*k4*E2T-a1*k1^2*d2*a3*d3*a4*d4-a1*k1^2*d2*a3*d3*a4*k4-a1*k1^2*k2*a3*d3*a4*d4-a1*k1^2*k2*a3*d3*a4*k4-a1*k1*a2*k2*a3*d3*a4*d4*E2T-a1*k1*a2*k2*a3*d3*a4*k4*E2T-a1*k1*a2*a3*d3*a4*d4*k4*E2T-a1*k1*a2*a3*d3*a4*k4^2*E2T-d1^2*a2*k2*a3*d3*a4^2*E2T+d1^2*d2*a3^2*a4^2*k4*WT+d1^2*k2*a3^2*a4^2*k4*WT-2*d1*k1*a2*k2*a3*d3*a4^2*E2T+2*d1*k1*d2*a3^2*a4^2*k4*WT+2*d1*k1*k2*a3^2*a4^2*k4*WT-k1^2*a2*k2*a3*d3*a4^2*E2T+k1^2*d2*a3^2*a4^2*k4*WT+k1^2*k2*a3^2*a4^2*k4*WT);
% *C4*Ws+
Coefficients(10) = (-a1^2*k1*a2*d3^2*d4^2-2*a1^2*k1*a2*d3^2*d4*k4-a1^2*k1*a2*d3^2*k4^2+a1*d1*k1*d2*a3*d3*a4*d4+a1*d1*k1*d2*a3*d3*a4*k4+a1*d1*k1*k2*a3*d3*a4*d4+a1*d1*k1*k2*a3*d3*a4*k4+a1*d1*a2*d3^2*a4*d4*k4+a1*d1*a2*d3^2*a4*k4^2-a1*d1*d2*a3*d3*a4*d4*k4-a1*d1*d2*a3*d3*a4*k4^2-a1*d1*k2*a3*d3*a4*d4*k4-a1*d1*k2*a3*d3*a4*k4^2+a1*k1^2*d2*a3*d3*a4*d4+a1*k1^2*d2*a3*d3*a4*k4+a1*k1^2*k2*a3*d3*a4*d4+a1*k1^2*k2*a3*d3*a4*k4+a1*k1*a2*d3^2*a4*d4*k4+a1*k1*a2*d3^2*a4*k4^2-a1*k1*d2*a3*d3*a4*d4*k4-a1*k1*d2*a3*d3*a4*k4^2-a1*k1*k2*a3*d3*a4*d4*k4-a1*k1*k2*a3*d3*a4*k4^2);
% *C4*E1T+
Coefficients(11) = (-a1^2*k1*a2*d3^2*d4^2*E2T-2*a1^2*k1*a2*d3^2*d4*k4*E2T-a1^2*k1*a2*d3^2*k4^2*E2T+a1*d1*k1*d2*a3*d3*a4*d4*WT+a1*d1*k1*d2*a3*d3*a4*k4*WT+a1*d1*k1*k2*a3*d3*a4*d4*WT+a1*d1*k1*k2*a3*d3*a4*k4*WT-a1*d1*a2*k2*d3^2*a4*d4*E2T-a1*d1*a2*k2*d3^2*a4*k4*E2T+a1*k1^2*d2*a3*d3*a4*d4*WT+a1*k1^2*d2*a3*d3*a4*k4*WT+a1*k1^2*k2*a3*d3*a4*d4*WT+a1*k1^2*k2*a3*d3*a4*k4*WT-a1*k1*a2*k2*d3^2*a4*d4*E2T-a1*k1*a2*k2*d3^2*a4*k4*E2T);
% *C4+
Coefficients(12) = (d1^2*a2*a3^2*a4^2*k4+2*d1*k1*a2*a3^2*a4^2*k4+k1^2*a2*a3^2*a4^2*k4);
% *Ws^3*E1T+
Coefficients(13) = (d1^2*a2*a3^2*a4^2*k4+2*d1*k1*a2*a3^2*a4^2*k4+k1^2*a2*a3^2*a4^2*k4);
% *Ws^2*E1T^2+
Coefficients(14) = (a1*d1*k1*a2*a3*d3*a4*d4+a1*d1*k1*a2*a3*d3*a4*k4+a1*k1^2*a2*a3*d3*a4*d4+a1*k1^2*a2*a3*d3*a4*k4-d1^2*a2*a3^2*a4^2*k4*WT+d1^2*a2*a3^2*a4^2*k4*E2T+d1^2*d2*a3^2*a4^2*k4+d1^2*k2*a3^2*a4^2*k4-2*d1*k1*a2*a3^2*a4^2*k4*WT+2*d1*k1*a2*a3^2*a4^2*k4*E2T+2*d1*k1*d2*a3^2*a4^2*k4+2*d1*k1*k2*a3^2*a4^2*k4-k1^2*a2*a3^2*a4^2*k4*WT+k1^2*a2*a3^2*a4^2*k4*E2T+k1^2*d2*a3^2*a4^2*k4+k1^2*k2*a3^2*a4^2*k4);
% *Ws^2*E1T+
Coefficients(15) = (-a1*d1*a2*k2*a3*d3*a4*d4*E2T-a1*d1*a2*k2*a3*d3*a4*k4*E2T-a1*k1*a2*k2*a3*d3*a4*d4*E2T-a1*k1*a2*k2*a3*d3*a4*k4*E2T+d1^2*a2*k2*a3^2*a4*d4*E2T+d1^2*a2*k2*a3^2*a4*k4*E2T+2*d1*k1*a2*k2*a3^2*a4*d4*E2T+2*d1*k1*a2*k2*a3^2*a4*k4*E2T+k1^2*a2*k2*a3^2*a4*d4*E2T+k1^2*a2*k2*a3^2*a4*k4*E2T);
% *Ws^2+
Coefficients(16) = (a1*d1*k1*a2*a3*d3*a4*d4+a1*d1*k1*a2*a3*d3*a4*k4+a1*k1^2*a2*a3*d3*a4*d4+a1*k1^2*a2*a3*d3*a4*k4-d1^2*a2*a3*d3*a4^2*k4+d1^2*d2*a3^2*a4^2*k4+d1^2*k2*a3^2*a4^2*k4-2*d1*k1*a2*a3*d3*a4^2*k4+2*d1*k1*d2*a3^2*a4^2*k4+2*d1*k1*k2*a3^2*a4^2*k4-k1^2*a2*a3*d3*a4^2*k4+k1^2*d2*a3^2*a4^2*k4+k1^2*k2*a3^2*a4^2*k4);
% *Ws*E1T^2+
Coefficients(17) = (-a1*d1*k1*a2*a3*d3*a4*d4*WT+a1*d1*k1*a2*a3*d3*a4*d4*E2T-a1*d1*k1*a2*a3*d3*a4*k4*WT+a1*d1*k1*a2*a3*d3*a4*k4*E2T+a1*d1*k1*d2*a3*d3*a4*d4+a1*d1*k1*d2*a3*d3*a4*k4+a1*d1*k1*k2*a3*d3*a4*d4+a1*d1*k1*k2*a3*d3*a4*k4-a1*k1^2*a2*a3*d3*a4*d4*WT+a1*k1^2*a2*a3*d3*a4*d4*E2T-a1*k1^2*a2*a3*d3*a4*k4*WT+a1*k1^2*a2*a3*d3*a4*k4*E2T+a1*k1^2*d2*a3*d3*a4*d4+a1*k1^2*d2*a3*d3*a4*k4+a1*k1^2*k2*a3*d3*a4*d4+a1*k1^2*k2*a3*d3*a4*k4+d1^2*a2*k2*a3*d3*a4^2*E2T-d1^2*d2*a3^2*a4^2*k4*WT-d1^2*k2*a3^2*a4^2*k4*WT+2*d1*k1*a2*k2*a3*d3*a4^2*E2T-2*d1*k1*d2*a3^2*a4^2*k4*WT-2*d1*k1*k2*a3^2*a4^2*k4*WT+k1^2*a2*k2*a3*d3*a4^2*E2T-k1^2*d2*a3^2*a4^2*k4*WT-k1^2*k2*a3^2*a4^2*k4*WT);
% *Ws*E1T+
Coefficients(18) = (a1*d1*a2*k2*a3*d3*a4*d4*WT*E2T+a1*d1*a2*k2*a3*d3*a4*k4*WT*E2T+a1*k1*a2*k2*a3*d3*a4*d4*WT*E2T+a1*k1*a2*k2*a3*d3*a4*k4*WT*E2T+d1^2*a2*k2*a3*d3*a4*d4*E2T+d1^2*a2*k2*a3*d3*a4*k4*E2T+2*d1*k1*a2*k2*a3*d3*a4*d4*E2T+2*d1*k1*a2*k2*a3*d3*a4*k4*E2T+k1^2*a2*k2*a3*d3*a4*d4*E2T+k1^2*a2*k2*a3*d3*a4*k4*E2T);
% *Ws+
Coefficients(19) = (-a1*d1*k1*d2*a3*d3*a4*d4*WT-a1*d1*k1*d2*a3*d3*a4*k4*WT-a1*d1*k1*k2*a3*d3*a4*d4*WT-a1*d1*k1*k2*a3*d3*a4*k4*WT-a1*k1^2*d2*a3*d3*a4*d4*WT-a1*k1^2*d2*a3*d3*a4*k4*WT-a1*k1^2*k2*a3*d3*a4*d4*WT-a1*k1^2*k2*a3*d3*a4*k4*WT);
% *E1T

% 1 C3 as a function of Ws, C4, E1T
Coefficients(20) = (d1*a3*a4+k1*a3*a4);
% *C3*Ws+
Coefficients(21) = (d1*d3*a4+k1*d3*a4);
% *C3+
Coefficients(22) = (d1*a3*a4+k1*a3*a4);
% *C4*Ws+
Coefficients(23) = (a1*d3*d4+a1*d3*k4);
% *C4+
Coefficients(24) = (-d1*a3*a4-k1*a3*a4);
% *Ws*E1T

% 2 C3 as a function of Ws,E1t,C4
Coefficients(25) = (d1*k2*d3*a4+k1*k2*d3*a4);
% *C3^2+
Coefficients(26) = (-d1*k2*a3*a4-d1*a3*a4*k4-k1*k2*a3*a4-k1*a3*a4*k4);
% *C3*C4*Ws+
Coefficients(27) = (-a1*k1*d3*d4-a1*k1*d3*k4);
% *C3*C4+
Coefficients(28) = (-d1*k2*a3*a4-k1*k2*a3*a4);
% *C3*Ws^2+
Coefficients(29) = (-d1*k2*a3*a4-k1*k2*a3*a4);
% *C3*Ws*E1T+
Coefficients(30) = (d1*k2*a3*a4*WT+k1*k2*a3*a4*WT);
% *C3*Ws+
Coefficients(31) = (-d1*k2*a3*d4-d1*k2*a3*k4-k1*k2*a3*d4-k1*k2*a3*k4);
% *C4*Ws

Coefficients(1:19) = Coefficients(1:19)./(a1^2*a2^2*a3^3*a4^2);

end

