function Y = density_Y(Y_T, Y_0, mu, T, sigma)
top = -(log(Y_T) - log(Y_0) - (mu - (1/2)*sigma.^2)*T).^2./(2*sigma.^2*T);
Y = 1./(Y_T .* sigma .* sqrt(2*pi*T)) .* exp(top);
end



