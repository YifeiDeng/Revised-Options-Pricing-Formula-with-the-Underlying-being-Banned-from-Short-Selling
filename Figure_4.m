% Replication of Figure #4

%% Initialized some variables
T = 1;
N = 10000;
dt = T/N;

sigma = 0.2;
mu = 0.1;
rho = 0.8;
space = 100;
Y_0 = linspace(5,20,space);
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
F_p = Forward(Y_0, K, mu, T);
delta = 1./(1-rho.^2);

% define matrix for plotting
C = zeros(1,space);
P = zeros(1,space);
F = zeros(1,space);
C_equ_risk= zeros(1,space);
P_equ_risk= zeros(1,space);
F_equ_risk = zeros(1,space);

for j = 1:length(Y_0)
    C(j) = delta/2 * log(integral(@(Y_T) density_Y(Y_T, Y_0(j), mu, T, sigma).*...
        exp(Payoff_fun(Y_T, K, mu, T, 'CALL')*(-1/delta)), 0, inf));
    P(j) = delta/2 * log(integral(@(Y_T) density_Y(Y_T, Y_0(j), mu, T, sigma).*...
        exp(Payoff_fun(Y_T, K, mu, T, 'PUT')*(1/delta)), 0, inf));
    F(j) = delta/2 * log(integral(@(Y_T) density_Y(Y_T, Y_0(j), mu, T, sigma).*...
        exp(exp(-mu*T) * (Y_T - K) * (-1/delta)), 0, inf));
    C_equ_risk(j) = 1/2 * C_bs(j) - C(j);
    P_equ_risk(j) = 1/2 * P_bs(j) + P(j);
    F_equ_risk(j) = 1/2 * F_p(j) - F(j);
end

P_C_parity = C_equ_risk - P_equ_risk - F_equ_risk;

hold on
plot(Y_0,P_C_parity,'-b')
hold off
grid on
xlabel 'Underlying Price'
ylabel 'Deviation of put-call parity'
title('Deviation of put-call parity against the Underlying Price')
