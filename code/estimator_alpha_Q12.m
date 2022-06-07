function alpha_n = estimator_alpha_Q12(Delta, X)
alpha_n = sum(X.^2) / sum(X(1:end-1) .* X(2:end));
alpha_n = 1 / Delta * log(alpha_n);
end

