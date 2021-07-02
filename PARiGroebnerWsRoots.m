function Coefficients = PARiGroebnerWsRoots(E2T,WT,k)
% this function calculates the coefficient of the groebner basis for Ws

a1 = k(1); d1 = k(2); k1 = k(3);
a2 = k(4); d2 = k(5); k2 = k(6);
a3 = k(7); d3 = k(8);
a4 = k(9); d4 = k(10); k4 = k(11);

K1 = (d1+k1)/a1; K2 = (d2+k2)/a2; Kd = d3/a3; K4 = (d4+k4)/a4; 

Coefficients = zeros(1,18); 
% % 6 values of Ws as a function of E1T
Coefficients(1) = k4^2*K1^2;
% % Ws^6*E1T
Coefficients(2) = -k2*k4*K1^2*E2T;
% % Ws^6
Coefficients(3) = k4^2*K1^2;
% % Ws^5*E1T^2
Coefficients(4) = k4*K1*(k1*Kd*K4*2 - k4*WT*K1 + k4*E2T*K1 + k4*Kd*K1 + 2*k4*K2*K1);
% % Ws^5*E1T
Coefficients(5) = k2*E2T*K1*(-k1*Kd*K4 - k4*K4*Kd - k2*E2T*K1 + k4*WT*K1 - k4*E2T*K1 + k4*K1*K4 - k4*K1*Kd - k4*K1*K2);
% % *Ws^5+
Coefficients(6) = 2*k4*K1*(k1*K4*Kd + k4*K1*K2);
% % *Ws^4*E1T^2+
Coefficients(7) = (k1^2*Kd^2*K4^2 - 2*k1*k4*WT*Kd*K4*K1 + 2*k1*k4*E2T*Kd*K4*K1 + 2*k1*k4*Kd^2*K4*K1 + 4*k1*k4*Kd*K4*K2*K1 - k2*k4*E2T*Kd*K4*K1 + 2*k2*k4*E2T*Kd*K1^2 - k4^2*WT*Kd*K1^2 + k4^2*E2T*Kd*K1^2 - 2*k4^2*WT*K1^2*K2 + k4^2*E2T*K1^2*K2 + 2*k4^2*Kd*K2*K1^2 + k4^2*K2^2*K1^2);
% % *Ws^4*E1T+
Coefficients(8) = (-k1*k2*E2T*Kd^2*K4^2 + k1*k2*WT*E2T*Kd*K4*K1 - k1*k2*E2T^2*Kd*K4*K1 + k1*k2*E2T*Kd*K4^2*K1 - k1*k2*E2T*Kd^2*K4*K1 - k1*k2*E2T*Kd*K4*K2*K1 - k2^2*E2T^2*Kd*K4*K1 + k2*k4*WT*E2T*Kd*K4*K1...
    - k2*k4*E2T^2*Kd*K4*K1 - k2*k4*E2T*Kd^2*K4*K1 - k2*k4*E2T*Kd*K4*K2*K1 + 2*k2*k4*E2T*Kd*K4*K1^2 + k2*k4*E2T*K4*K2*K1^2 + k2*k4*WT*E2T*K2*K1^2 - k2*k4*E2T*Kd*K2*K1^2 - 2*k2^2*E2T^2*Kd*K1^2 + k2*k4*WT*E2T*Kd*K1^2 - k2*k4*E2T^2*Kd*K1^2);
% % *Ws^4+
% k4 term
Coefficients(9) = (k4^2*K1^2*K2^2 + k1^2*Kd^2*K4^2 + 4*k1*k4*Kd*K1*K2*K4);
% % *Ws^3*E1T^2+
Coefficients(10) = (-k1^2*WT*Kd^2*K4^2 + k1^2*E2T*Kd^2*K4^2 + k1^2*Kd^3*K4^2 + 2*k1^2*Kd^2*K4^2*K2 - k1*k2*E2T*Kd^2*K4^2 + 2*k1*k2*E2T*Kd^2*K4*K1 - 2*k1*k4*WT*Kd^2*K4*K1 + 2*k1*k4*E2T*Kd^2*K4*K1 - 4*k1*k4*WT*Kd*K4*K2*K1 + 2*k1*k4*E2T*Kd*K4*K2*K1 + 4*k1*k4*Kd^2*K4*K2*K1 + 2*k1*k4*Kd*K4*K2^2*K1 + k2*k4*E2T*Kd^2*K4*K1 - k2*k4*E2T*Kd*K4*K2*K1 - 2*k4^2*WT*Kd*K2*K1^2 + 2*k2*k4*E2T*Kd*K2*K1^2 + k4^2*E2T*Kd*K2*K1^2 - k4^2*WT*K2^2*K1^2 + k4^2*Kd*K2^2*K1^2);
% % *Ws^3*E1T+
Coefficients(11) = Kd*(k1*k2*WT*E2T*Kd*K4^2 - k1*k2*E2T^2*Kd*K4^2 - k1*k2*E2T*Kd^2*K4^2 - k1*k2*E2T*Kd*K2*K4^2 + k1*k2*WT*E2T*Kd*K4*K1 - k1*k2*E2T^2*Kd*K4*K1 + 2*k1*k2*E2T*Kd*K4^2*K1 + k1*k2*WT*E2T*K4*K2*K1 + k1*k2*E2T*K4^2*K2*K1 - k1*k2*E2T*Kd*K4*K2*K1 ...
    - 3*k2*k2*E2T^2*Kd*K4*K1 + k2*k4*WT*E2T*Kd*K4*K1 - k2*k4*E2T^2*Kd*K4*K1 + k2*k4*WT*E2T*K4*K2*K1 - k2*k4*E2T*Kd*K4*K2*K1 + k2*k4*E2T*Kd*K4*K1^2 + k2*k4*WT*E2T*K2*K1^2 + 2*k2*k4*E2T*K4*K2*K1^2);
% % *Ws^3+
Coefficients(12) = 2*k1*Kd*K2*K4*(k1*Kd*K4 + k4*K2*K1);
% % *Ws^2*E1T^2+
Coefficients(13) = Kd*(-k1^2*WT*Kd^2*K4^2 + k1^2*E2T*Kd^2*K4^2 - 2*k1^2*WT*Kd*K4^2*K2 + k1^2*E2T*Kd*K4^2*K2 + 2*k1^2*Kd^2*K4^2*K2 + k1^2*Kd*K4^2*K2^2 + k1*k2*E2T*Kd^2*K4^2 - k1*k2*E2T*Kd*K4^2*K2 + 2*k1*k2*E2T*Kd*K1*K2*K4 - 4*k1*k4*WT*Kd*K4*K2*K1 + 2*k1*k4*E2T*K1*K2*Kd*K4 - 2*k1*k4*WT*K4*K2^2*K1 + 2*k1*k4*Kd*K4*K2^2*K1 + k2*k4*E2T*Kd*K4*K2*K1 - k4^2*WT*K2^2*K1^2);
% % *Ws^2*E1T+
Coefficients(14) = k2*K4*Kd^2*(k1*WT*E2T*Kd*K4 - k1*E2T^2*Kd*K4 + k1*WT*E2T*K4*K2 - k1*E2T*K2*Kd*K4 - k2*E2T^2*Kd*K4 + k1*E2T*Kd*K4*K1 + k1*WT*E2T*K2*K1 + 2*k1*E2T*K4*K2*K1 + k4*WT*E2T*K2*K1 + k4*E2T*K1^2*K2);
% % *Ws^2+
Coefficients(15) = k1^2*K2^2*K4^2*Kd^2;
% % *Ws*E1T^2+
Coefficients(16) = k1*K2*K4*Kd^2*(-2*k1*Kd*WT*K4 + k1*Kd*E2T*K4 - k1*WT*K2*K4 + k1*K2*K4*Kd + k2*E2T*Kd*K4 - 2*k4*WT*K1*K2);
% % *Ws*E1T+
Coefficients(17) = k1*k2*E2T*K2*Kd^3*K4^2*(WT + K1);
% % *Ws+ 
Coefficients(18) = -k1^2*WT*Kd^3*K2^2*K4^2;
% % *E1T


end
 