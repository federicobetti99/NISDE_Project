function g = G_2D_Q15(x, X, delta)
    g = sum((...
    X(2:end).^2 + X(2:end) ...
    - exp(-2 * x(1) * delta) * X(1:end-1).^2 ...
    - exp(- x(1) * delta) * X(1:end-1) ...
    - (1 - exp(-2 * x(1) * delta)) * x(2) / x(1)) ...
    .* [X(1:end-1).^2; X(1:end-1)], 2);
end