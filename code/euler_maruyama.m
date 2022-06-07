function X = euler_maruyama(x0, alpha, sigma, delta_t, T)
    N = T/delta_t;
    X = zeros(1, N+1);
    X(1) = x0;
    for n = 1:N
        X(n+1) = X(n) - alpha * X(n) * delta_t + sqrt(2*sigma) * sqrt(delta_t) * randn(1);
    end
end

