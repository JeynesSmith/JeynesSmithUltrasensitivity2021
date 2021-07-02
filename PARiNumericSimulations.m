% This function is used to simulate the PAR-i system, and find steady
% states across a domain of input values defind by x. The initial
% conditions are defined by init which is dependent on x, Wtot, C1, C2, C3,
% C4, W, Ws, and E2. The parameters are defined by a1, d1, ..., d4, k4.
% Simultions are then combined in a single figure to generate figures in
% paper

% Set the initial values
Wtot = 100; 
C1 = 0; C2 = 0; C3 = 0; C4 = 0;
W = Wtot; Ws = 0;
E2 = 20;

% define input 
x= 0:1:40;

% Set the model parameters
a1 = 1; d1 = 1; k1 = 1;
a2 = 1; d2 = 1; k2 = 1;
a3 = 0.1; d3 = 10;
a4 = 1; d4 = 1; k4 = 1;
k = [a1,d1,k1,a2,d2,k2,a3,d3,a4,d4,k4];

% initialise vectors
Wplot = zeros(size(x));
Wsplot = zeros(size(x));

% Start Lopping Process
tic
for i = 1:length(x)
    
    % update initial conditions
    E1 = x(i);
    init = [W Ws E1 E2 C1 C2 C3 C4];
    
    % Set the time domain
    tspan = [0 1e6];
    
    % Perform the numerical integration
    [t,u] = ode23s(@(t,u) odesys(t,u,k), tspan, init);
    
    % Check Convergence
    if sum(abs(u(end,:) - u(end-1,:)) > 1e-2) >0
        warning('Timestep too small')
        disp((u(end,:) - u(end-1,:) > 1e-2) .* ( u(end,:) - u(end-1,:)))
        disp(E1)
    end
    
    % Create vectors for plots
    Wplot(i) = u(end,1)/Wtot;
    Wsplot(i) = u(end,2)/Wtot;
    
end
toc

% Plot the concentration of Ws across input values
figure
plot(x,Wsplot,'-','LineWidth',2,'color',[1,0.5,0])
ylim([0,1])
xlabel('Input (E_{1tot})','FontSize',14)
ylabel('Output (W*)','FontSize',14)

%% ODE System
function eqns = odesys(t,u,k)
eqns = zeros(8,1); % To start with we have six empty equations
% Using u = [W  Ws E1 E2 C1 C2 C3 C4] 
% Using k = [a1,d1,k1,a2,d2,k2,a3,d3,a4,d4,k4]
%           [ 1  2  3  4  5  6  7  8  9 10 11
eqns(1) = k(2)*u(5) + k(6)*u(6) + k(10)*u(8) - k(1)*u(1)*u(3) - k(9)*u(1)*u(7);
eqns(2) = k(5)*u(6) + k(3)*u(5) + k(8)*u(7) + k(11)*u(8) - k(4)*u(2)*u(4) - k(7)*u(2)*u(3);
eqns(3) = k(2)*u(5) + k(3)*u(5) + k(8)*u(7) - k(1)*u(1)*u(3) - k(7)*u(2)*u(3);
eqns(4) = (k(5) + k(6))*u(6) - k(4)*u(2)*u(4);
eqns(5) = k(1)*u(1)*u(3) - (k(2) + k(3))*u(5);
eqns(6) = k(4)*u(2)*u(4) - (k(5) + k(6))*u(6);
eqns(7) = k(7)*u(2)*u(3) + (k(10) + k(11))*u(8) - k(8)*u(7) - k(9)*u(1)*u(7);
eqns(8) = k(9)*u(1)*u(7) - (k(10) + k(11))*u(8);
end