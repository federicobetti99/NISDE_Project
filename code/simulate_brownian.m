function W = simulate_brownian(a, b, N)
    sqrt_delta = sqrt((b-a) / N);
    W = [0, cumsum(sqrt_delta * randn(1,N))];
end

