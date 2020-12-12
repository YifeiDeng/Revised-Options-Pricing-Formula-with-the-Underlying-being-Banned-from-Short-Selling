function Payoff = Payoff_fun(S_T, K, r, T, option_type)

switch option_type
    case 'CALL'
        % Payoff function for European Call Option
        Payoff = exp(-r * T) * max(S_T - K, 0);
    case 'PUT'
        % Payoff function for European Put Option
        Payoff = exp(-r * T) * max(K - S_T, 0);
end

