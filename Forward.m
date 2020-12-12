function F = Forward(S_0, K, r, T)
    F = S_0 - exp(-r * T) * K;
end