% Replication of Figure #1

%% Initialized some variables
T = 1;
N = 10000;
dt = T/N;

mu = 0.1;
sigma = 0.2;
rho = [0, 0.5];
space = 20;
Y_0 = linspace(5,15,space);
K = 10;

%% simulate Y_t

dW_1 = normrnd(0,1,[N,space])*sqrt(dt);

Y_t = zeros(N+1,space);
Y_t(1,:) = Y_0;

% Start simulating
for i = 2:(N+1)
    dY_t = Y_t(i-1,:) * mu * dt + sigma * dW_1(i-1,:);
    Y_t(i,:) = Y_t(i-1,:) + dY_t;
end

Y_T = Y_t(end,:);

%% simulate Call and Put Option price for European Option

% for rho not equal to -1 or +1

% Initialized some known values
C_bs = Exact_BS(Y_0, K, mu, sigma, T, 'CALL');
P_bs = Exact_BS(Y_0, K, mu, sigma, T, 'PUT');
delta = 1./(1-rho.^2);

% define matrix for plotting
C = zeros(space,3);
P = zeros(space,3);
C_equ_risk= zeros(space,3);
P_equ_risk= zeros(space,3);

for i = 1:length(Y_0)
    for j = 1:2
        C(i,j) = delta(j)/2 * log(integral(@(Y_T) density_Y(Y_T, Y_0(i), mu, T, sigma).*...
                                  exp(Payoff_fun(Y_T, K, mu, T, 'CALL')*(-1/delta(j))), 0, inf));
        P(i,j) = delta(j)/2 * log(integral(@(Y_T) density_Y(Y_T, Y_0(i), mu, T, sigma).*...
                                  exp(Payoff_fun(Y_T, K, mu, T, 'PUT')*(1/delta(j))), 0, inf));
        C_equ_risk(i,j) = 1/2 * C_bs(i) - C(i,j);
        P_equ_risk(i,j) = 1/2 * P_bs(i) + P(i,j);
    end
end

C_equ_risk(:,3) = C_bs;
P_equ_risk(:,3) = P_bs;

hold on
plot(Y_0, C_equ_risk(:,1),'-.b')
plot(Y_0, C_equ_risk(:,2),'--b')
plot(Y_0, C_equ_risk(:,3),':b')
plot(Y_0, P_equ_risk(:,1),'-or')
plot(Y_0, P_equ_risk(:,2),'-+r')
plot(Y_0, P_equ_risk(:,3),'-*r')
hold off
legend show
legend('Call with \rho = 0 (Guo-Zhu)', 'Call with \rho = 0.5',...
       'Call with \rho = 1 (Black-Scholes)', 'Put with \rho = 0 (Guo-Zhu)',...
       'Put with \rho = 0.5', 'Put with \rho = 1 (Black-Scholes)')
xlabel 'Underlying Price'
ylabel 'Option Price'
legend('Location','north')
title('European Call and Put prices against the Underlying Price')



