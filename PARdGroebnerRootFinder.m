% This function is used to calculate every steady-state across a domain of
% input values defind by iter in the PAR-d system. The total abundances are
% defined by Wtot and E2T. The parameters are defined by a1, d1, ..., d3, k3.
% Roots are calculated using a Groebner basis, where the coefficients are
% given by PARiGroebnerWsRoots.m, and PARiGroebnerWRoots.m. Stability of 
% the roots is then determined by using a linear stability test. All roots 
% are then plotted in a figure where the filling of the marker is 
% determined by the stability: solid = stable, empty = unstable. 

% Set the total abundances
Wtot = 100;
E2T = 20;

% Set the model parameters
a1 = 1; d1 = 1; k1 = 1; 
a2 = 1; d2 = 1; k2 = 1;
a3 = 0.01; d3 = 1; k3 = 1;
k = [a1,d1,k1,a2,d2,k2,a3,d3,k3];

% initialise roots
EveryRoot = [];
for iter = 0:1:40 % loop over input values
    E1T = iter;
    
    % Find the values of Ws
    A = PARdGroebnerWsCoefs(k,Wtot,E2T);
    WsRoots = roots([A(1),A(2)*E1T+A(3),A(4)*E1T+A(5),A(6)*E1T+A(7),A(8)*E1T+A(9),A(10)*E1T]);
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
    
    % find the values of W
    A = PARdGroebnerWCoefs(k,Wtot,E2T);
    Combo(:,3) = (A(2).*Combo(:,2).^4 + A(3).*Combo(:,1).*Combo(:,2).^3 + A(4).*Combo(:,2).^3 + A(5).*Combo(:,1).*Combo(:,2).^2 + A(6).*Combo(:,2).^2 + A(7).*Combo(:,1).*Combo(:,2) + A(8).*Combo(:,2) + A(9).*Combo(:,1) + A(10))./(-A(1));
    % remove negative roots
    Combo = Combo(Combo(:,3)>=0,:);
    
    % find the values of C3    
    Combo(:,4) = a3./(d3+k3).*Combo(:,2).*Combo(:,3);
    
    % find the values of C2
    Combo(:,5) = ((2*k1-k3).*Combo(:,4) + k1.*Combo(:,3) + k1.*Combo(:,2) - k1*Wtot)./(-k1-k2);
    % remove negative roots
    Combo = Combo(Combo(:,5)>=0,:);
    
    % find the values of C1, E1, E2
    Combo(:,6:7) = [Wtot-Combo(:,2)-Combo(:,3)-2.*Combo(:,4)-Combo(:,5),E2T-Combo(:,5)];
    Combo(:,8) = Combo(:,1)-Combo(:,6);
    Combo = Combo(all(Combo(:,6:8)>=0,2),:);
    
    % store all valid roots
    EveryRoot = [EveryRoot; Combo];
    
end

% plot roots
figure(1), clf, hold on
stabPoints = []; unstabPoints = [];
for i = 1:size(EveryRoot,1)
    % check stability of roots
    % EveryRoot = [E1T Ws W C3 C2 C1 E2 E1]
    J = [
        -a1*EveryRoot(i,8)-a3*EveryRoot(i,2), -a3*EveryRoot(i,3), -a1*EveryRoot(i,3), 0, d1, k2, d3
        -a3*EveryRoot(i,2), -a2*EveryRoot(i,7)-a3*EveryRoot(i,3), 0, -a2*EveryRoot(i,2), k1, d2, d3+2*k3
        -a1*EveryRoot(i,8), 0, -a1*EveryRoot(i,3), 0, d1+k1, 0, 0
        0, -a2*EveryRoot(i,7), 0, -a2*EveryRoot(i,2), 0, d2+k2, 0
        a1*EveryRoot(i,8), 0, a1*EveryRoot(i,3), 0, -d1-k1, 0, 0
        0, a2*EveryRoot(i,7), 0, a2*EveryRoot(i,2), 0, -d2-k2, 0
        a3*EveryRoot(i,2), a3*EveryRoot(i,3), 0, 0 ,0, 0, -d3-k3
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

