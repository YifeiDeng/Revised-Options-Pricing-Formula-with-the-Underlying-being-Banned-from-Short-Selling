% Replication of Figure #2

%% Initialized some variables
T = 1;
N = 10000;
dt = T/N;

sigma = 0.2;
mu = 0.1;
space = 19;
rho_1 = linspace(0,1,space+1);
rho = rho_1(1:(end-1));
Y_0 = 10;
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
C = zeros(1,space);
P = zeros(1,space);
C_equ_risk= zeros(1,space);
P_equ_risk= zeros(1,space);

for j = 1:length(rho)
    C(j) = delta(j)/2 * log(integral(@(Y_T) density_Y(Y_T, Y_0, mu, T, sigma).*...
        exp(Payoff_fun(Y_T, K, mu, T, 'CALL')*(-1/delta(j))), 0, inf));
    P(j) = delta(j)/2 * log(integral(@(Y_T) density_Y(Y_T, Y_0, mu, T, sigma).*...
        exp(Payoff_fun(Y_T, K, mu, T, 'PUT')*(1/delta(j))), 0, inf));
    C_equ_risk(j) = 1/2 * C_bs - C(j);
    P_equ_risk(j) = 1/2 * P_bs + P(j);
end

C_equ_risk = [C_equ_risk, C_bs];
P_equ_risk = [P_equ_risk, P_bs];

hold on
plot(rho_1,C_equ_risk,'--b')
plot(rho_1,P_equ_risk,'-or')
hold off
legend('Call', 'Put')
legend show
grid on
legend('Location','north')
xlabel 'Correlation'
ylabel 'Option Price'
title('European Call and Put Prices against the ban-dillution factor')



