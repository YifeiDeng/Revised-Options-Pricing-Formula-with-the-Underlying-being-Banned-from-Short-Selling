% Exact solution for a European type option in Black-Scholes
% equations
%--------------------------------------------------------------------------
% INPUTS:
%
%   K:      strike price
%   S_0:    the spot stock price
%   r:      risk-free interest rate
%   T:      time to maturity
%   sigma:  volatility
%   option_type:  option types, type either'CALL' or 'PUT'
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% V: option value
%--------------------------------------------------------------------------
function V = Exact_BS(S_0, K, r, sigma, T, option_type)
% Determine d1 and d2
d1 = (log(S_0./K) + (r + sigma^2/2)*T)/(sigma * sqrt(T));
d2 = d1 - sigma * sqrt(T);
switch option_type 
    case 'CALL'
        % Exact European Call Option solution using BS metron model
        V = S_0 .* normcdf(d1) - K * exp(-r * T) * normcdf(d2);
    case 'PUT'
        % Exact European Put Option solution using BS metron model
        V = K * exp(-r * T) * normcdf(-d2) - S_0 .* normcdf(- d1);
    
end