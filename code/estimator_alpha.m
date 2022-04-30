function alpha_n = estimator_alpha(Delta, X)
alpha_n = sum(X(1:end-1) .* X(2:end)) / sum(X.^2);
alpha_n = 1/(Delta) * (1 - alpha_n);
end

