%% parameters
T = 1e3; % final time
M = 1e4; % averages
h = 0.01; % step size for EM
Deltas = 2.^(-(0:7)); 
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition

%% compute estimators of \alpha and \sigma
error_estimates_alpha = zeros(length(Deltas), 1);
error_estimates_sigma = zeros(length(Deltas), 1);
for i=1:length(Deltas)
    W = simulate_brownian(0, T, T/h-1);
    X = euler_maruyama(x0, alpha, sigma, h, W);
    X_observed = X(1:(Deltas(i)/h):end);
    error_estimates_alpha(i) = abs(estimator_alpha(Deltas(i), X_observed)-alpha);
    error_estimates_sigma(i) = abs(estimator_sigma(Deltas(i), X_observed)-sigma);
end

%% plot of estimated \alpha and \sigma against the true ones
figure()
semilogx(Deltas, error_estimates_alpha, 'Linewidth', 1.0);
hold on
semilogx(Deltas, error_estimates_sigma, 'Linewidth', 1.0);
legend("$\vert \hat{\alpha_N}^{\Delta} - \alpha \vert$", "$ \vert \hat{\sigma_N}^{\Delta} - \sigma \vert$", ...
       "interpreter", "latex", "Location", "Best", "FontSize", 10);
title("Convergence of the estimators $\hat{\alpha_N}^{\Delta}$ and $\hat{\sigma_N}^{\Delta}$",  "interpreter", "latex");
saveas(gcf, "plot_Q7", "epsc");

