function sigma_n = estimator_sigma(Delta, X)
sigma_n = sum((X(2:end)-X(1:end-1)).^2);
sigma_n = 1/(2*Delta*length(X)) * sigma_n;
end

