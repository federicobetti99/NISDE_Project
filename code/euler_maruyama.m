function X = euler_maruyama(x0, alpha, sigma, delta_t, W)
    N = length(W) - 1;
    X = zeros(1, N+1);
    X(1) = x0;
    for n = 1:N
        X(n+1) = X(n) - alpha * X(n) * delta_t + sqrt(2*sigma) * (W(n+1) - W(n));
    end
end

