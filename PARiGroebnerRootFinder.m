% This function is used to calculate every steady-state across a domain of
% input values defind by iter in the PAR-i system. The total abundances are
% defined by Wtot and E2T. The parameters are defined by a1, d1, ..., d4, k4.
% Roots are calculated using a Groebner basis, where the coefficients are
% given by PARiGroebnerWsRoots.m, PARiGroebnerC4Roots.m, and
% PARiGroebnerC3Roots.m. Stability of the roots is then determined by using
% a linear stability test. All roots are then plotted in a figure where the
% filling of the marker is determined by the stability: solid = stable,
% empty = unstable. 

% Set the total abundances
Wtot = 100;
E2T = 20;

% Set the model parameters
a1 = 1; d1 = 1; k1 = 1;
a2 = 1; d2 = 1; k2 = 1;
a3 = 0.1; d3 = 10;
a4 = 1; d4 = 1; k4 = 10;
k = [a1,d1,k1,a2,d2,k2,a3,d3,a4,d4,k4];

% initialise roots
EveryRoot = [];
for iter = 0:40 % loop over input values
    E1T = iter;
    
    % Find the values of Ws
    A = PARiGroebnerWsRoots(E2T,Wtot,k);
    WsRoots = roots([A(1)*E1T + A(2), A(3)*E1T^2 + A(4)*E1T + A(5), A(6)*E1T^2 + A(7)*E1T + A(8), A(9)*E1T^2 + A(10)*E1T + A(11), A(12)*E1T^2 + A(13)*E1T + A(14), A(15)*E1T^2 + A(16)*E1T + A(17), A(18)*E1T]);
    % remove complex and negative roots
    Keep = zeros(size(WsRoots));
    for i = 1:length(Keep)
        if imag(WsRoots(i)) == 0 && real(WsRoots(i))>=0
            Keep(i) = 1;
        end
    end
    WsRoots = WsRoots(logical(Keep));
    % Store values
    Combo = [ones(size(WsRoots)).*E1T,WsRoots,zeros(size(WsRoots))];
    
    % Find the values of C4
    A = PARiGroebnerC4Roots(E2T,Wtot,k);
    Combo(:,3) = (Combo(:,2).^5.*(Combo(:,1).*A(2) + A(3)) + Combo(:,2).^4.*(Combo(:,1).^2.*A(4) + Combo(:,1).*A(5) + A(6)) + Combo(:,2).^3.*(Combo(:,1).^2.*A(7) + Combo(:,1).*A(8) + A(9))+ Combo(:,2).^2.*(Combo(:,1).^2.*A(10) + Combo(:,1).*A(11) + A(12))+ Combo(:,2).*(Combo(:,1).^2.*A(13) + Combo(:,1).*A(14) + A(15)) + A(16).*Combo(:,1))./(-A(1));
    % Store values
    Combo = Combo(Combo(:,3)>=0,:);
    
    % Find the values of C3
    A = PARiGroebnerC3Roots(E2T,Wtot,k);
    C3RootsEq1 = (Combo(:,3).^2.*(Combo(:,2).^2.*A(3) + Combo(:,2).*A(4) + A(5)) + Combo(:,3).*(Combo(:,2).^3.*A(6) + Combo(:,2).^2.*A(7) + Combo(:,2).*(Combo(:,1).*A(8)+A(9)) + Combo(:,1).*A(10) + A(11)) + Combo(:,2).^3.*Combo(:,1).*A(12) + Combo(:,2).^2.*(Combo(:,1).^2.*A(13) + Combo(:,1).*A(14) + A(15)) + Combo(:,2).*(Combo(:,1).^2.*A(16) + Combo(:,1).*A(17) + A(18)) + Combo(:,1).*A(19))./(-Combo(:,1).*A(1)-A(2));
    C3RootsEq2 = (Combo(:,3).*Combo(:,2)*A(22) + Combo(:,3).*A(23) + Combo(:,2).*Combo(:,1).*A(24))./(-Combo(:,2)*A(20) - A(21));
    ComboTempLength = size(Combo,1);
    Temp = zeros(ComboTempLength*4,9);
    for i = 1:ComboTempLength
        C3RootsEq3 = roots([A(25), Combo(i,3).*(Combo(i,2).*A(26) + A(27)) + Combo(i,2).^2.*A(28) + Combo(i,2).*(Combo(i,1).*A(29) + A(30)), Combo(i,3).*Combo(i,2).*A(31)]);
        Temp((i*4-3):(i*4),1:4) = [ones(4,3).*Combo(i,:) , [C3RootsEq1(i); C3RootsEq2(i); C3RootsEq3]];
    end
    % remove complex and negative roots
    Keep = zeros(size(Temp,1),1);
    for i = 1:length(Keep)
        if imag(Temp(i,4)) == 0 && real(Temp(i,4))>=0
            Keep(i) = 1;
        end
    end
    Combo = Temp(logical(Keep),:);
    
    % calculate the values of C2
    Combo(:,5) = ((-k1*a2*a3).*Combo(:,4).*Combo(:,2)+(-k1*a2*d3).*Combo(:,4)+(-k1*a2*a3+a2*a3*k4).*Combo(:,3).*Combo(:,2)+(k1*a2*a3).*Combo(:,2).*Combo(:,1)+(-a2*k2*a3*E2T).*Combo(:,2))./(-d2*k2*a3-k2^2*a3);
    Combo = Combo(Combo(:,5)>=0,:);
    % calculate the values of C1
    Combo(:,6) = ((-k2).*Combo(:,5)+(k4).*Combo(:,3))./(-k1);
    Combo = Combo(Combo(:,6)>=0,:);
    % calculate the values of E2, E1, W
    Combo(:,7:9) = [E2T-Combo(:,5), Combo(:,1)-Combo(:,6)-Combo(:,4)-Combo(:,3), Wtot-Combo(:,2)-Combo(:,6)-Combo(:,5)-Combo(:,4)-2.*Combo(:,3)];
    Combo = Combo(all(Combo(:,7:9)>=0,2),:);
    
    % store all valid roots
    EveryRoot = [EveryRoot; Combo];
    
end

% plot roots
figure(1), clf, hold on
stabPoints = []; unstabPoints = [];
for i = 1:size(EveryRoot,1)
    % check stability of roots
    % EveryRoot = [E1T Ws C4 C3 C2 C1 E2 E1 W]
    J = [
        -a1*EveryRoot(i,8)-a4*EveryRoot(i,4), 0, -a1*EveryRoot(i,9), 0, d1, k2, -a4*EveryRoot(i,9), d4;
        0, -a2*EveryRoot(i,7)-a3*EveryRoot(i,8), -a3*EveryRoot(i,2), -a2*EveryRoot(i,2), k1, d2, d3, k4;
        -a1*EveryRoot(i,8), -a3*EveryRoot(i,8), -a1*EveryRoot(i,9)-a3*EveryRoot(i,2), 0, d1+k1, 0, d3, 0;
        0, -a2*EveryRoot(i,7), 0, -a2*EveryRoot(i,2), 0, d2+k2, 0, 0;
        a1*EveryRoot(i,8), 0, a1*EveryRoot(i,9), 0, -d1-k1, 0, 0, 0;
        0, a2*EveryRoot(i,7), 0, a2*EveryRoot(i,2), 0, -d2-k2, 0, 0;
        -a4*EveryRoot(i,4), a3*EveryRoot(i,8), a3*EveryRoot(i,2), 0, 0, 0, -d3-a4*EveryRoot(i,9), d4+k4;
        a4*EveryRoot(i,4), 0, 0, 0 ,0, 0, a4*EveryRoot(i,9), -d4-k4
        ];
    TempEvals = eig(J);
    % plot and add to a vector based on stability. 
    if sum(real(TempEvals)>=1e-5) > 0
        scatter(EveryRoot(i,1),EveryRoot(i,2)./Wtot,'o','MarkerEdgeColor',[1 0.5 0],'LineWidth',2)
        unstabPoints = [unstabPoints; EveryRoot(i,1),EveryRoot(i,2)];
    else
        scatter(EveryRoot(i,1),EveryRoot(i,2)./Wtot,'o','MarkerEdgeColor',[1 0.5 0],'LineWidth',2,'MarkerFaceColor',[1 0.5 0])
        stabPoints = [stabPoints; EveryRoot(i,1),EveryRoot(i,2)];
    end
    
end
ax = gca;
ax.FontSize=14;
ylabel('Output (W*)','fontsize',17), xlabel('Input (E_{1 tot})','fontsize',17),ylim([0 1])

