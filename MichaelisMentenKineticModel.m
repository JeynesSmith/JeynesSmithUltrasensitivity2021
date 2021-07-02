% this function is used to examine the model based on Michaelis-Menten
% kinetics. The parameters K1, K2, k1, and k2 are first selected, along
% with the total abundances E2T and Wtot. The first figures examines the
% positive and negative contributions to the equation for a few input 
% values given by E1T. The second figure examines the steady states of the
% system across the same range of input values, but at a much higher
% sampling rate. We use three methods to test stability of the roots, which
% are toggled on and off by commenting them out. These methods have varying
% success depending on the parmeter set used. 

% define parameters
K1 = 1e-5;
K2 = 1e5;
k1 = 1;
k2 = 1;
%define total abundances
E2T = 20;
Wtot = 100;
% define input domain
% E1T = [0,2.5];
E1T = linspace(0,40/K2,6);
WsRange = 0:1:Wtot;

% figure 1
figure(1), clf
xlabel('Concentration of W*','fontsize',20)
ylabel('dW*/dt','fontsize',20)
hold on
% calculate negative contribution
NWs = k2.*E2T.*WsRange./(K2+WsRange);
Nmax = max(NWs);
% calculate positive contribution and plot
for i = 1:length(E1T)
    PWs = k1.*E1T(i).*WsRange.*(Wtot-WsRange)./(K1+Wtot-WsRange);
    plot(WsRange,PWs,'-','color',[1 0.5 0].*i./length(E1T),'linewidth',2)
    [M,ind] = max(PWs);
    Nmax = max([Nmax, M]);
end
plot(WsRange,NWs,'-','color',[0 0.5 1],'linewidth',3)
ylim(1.1.*[0, Nmax])

% figure 2
% increase sampling rate of input
E1Range = linspace(E1T(1),E1T(end),50);
x = [];
% calculate roots
for i = 1:length(E1Range)
    temproots = roots([-k1*E1Range(i),k1*E1Range(i)*(Wtot-K2)+k2*E2T,k1*E1Range(i)*Wtot*K2 - k2*E2T*(K1+Wtot),0]);
    x = [x; E1Range(i).*ones(size(temproots)),temproots];
end
% remove imaginary, negative or roots which are larger than the total
% substrate
Keep = zeros(length(x));
for i = 1:length(Keep)
    if imag(x(i,2)) == 0 && real(x(i,2))>=0 && real(x(i,2))<=Wtot
        Keep(i) = 1;
    end
end
x = x(logical(Keep),:);

% check stability and plot
tol = 10^-3; col = [1,0.5,0]; deviation = 1e-2;
deri = @(y,z) (-k1*z*y)/(K1+Wtot-y) + (-k1*z*y*(Wtot-y))/(K1+Wtot-y)^2 + (-k1*z*(Wtot-y))/(K1+Wtot-y) - k2*E2T/(K2+y) + k2*E2T*y/(K2+y)^2;
figure(2), clf, hold on
unstabPoints = []; stabPoints = [];
for i = 1:length(x)
    e = Jacobian(Wtot-x(i,2),x(i,2),k1,K1,k2,K2,x(i,1),E2T);
    if sum(real(e)>tol) > 0 % check 1
%     if deri(x(i,2),x(i,1)) > -tol % check 2
%     if func(x(i,2)-deviation,x(i,1),E2T,k1,k2,K1,K2,Wtot) < 0 ||
%     func(x(i,2)+deviation,x(i,1),E2T,k1,k2,K1,K2,Wtot) > 0 % check 3
        scatter(x(i,1),x(i,2)./100,'o','MarkerEdgeColor',[0,0.5,1],'LineWidth',2)
        unstabPoints = [unstabPoints; x(i,1),x(i,2)];
    else
        scatter(x(i,1),x(i,2)./100,'o','MarkerEdgeColor',col,'LineWidth',2,'MarkerFaceColor',col)
        stabPoints = [stabPoints; x(i,1),x(i,2)];
    end
end

% check 2 and 3
function eqns = func(W,E1T,E2T,k1,k2,K1,K2,Wtot)
eqns = k1.*E1T.*W.*(Wtot-W)./(K1+Wtot-W) - k2.*E2T.*W./(K2+W);
end
% check 1
function evalues = Jacobian(W,Ws,k1,K1,k2,K2,E1T,E2T)
J = [ -k1*E1T*Ws*K1/(K1+W)^2, k2*E2T*K2/(K2+Ws)^2 - k1*E1T*W/(K1+W);
    k1*E1T*Ws*K1/(K1+W)^2, -k2*E2T*K2/(K2+Ws)^2 + k1*E1T*W/(K1+W);];
evalues = eig(J);
end
