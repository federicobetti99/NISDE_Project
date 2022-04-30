%% parameters
T = 1e3; % final time
M = 1e4; % averages
h = 0.001; % step size for EM
Deltas = 2.^(-(0:7)); 
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition

%% estimate of \alpha using Q11)
error_estimates_alpha = zeros(length(Deltas), 1);
for i=1:length(Deltas)
    W = simulate_brownian(0, T, T/h-1);
    X = euler_maruyama(x0, alpha, sigma, h, W);
    X_observed = X(1:(Deltas(i)/h):end);
    error_estimates_alpha(i) = abs(estimator_alpha_Q11(Deltas(i), X_observed)-alpha);
end

%% plot of estimated \alpha
figure()
semilogx(Deltas, error_estimates_alpha, 'Linewidth', 1.0);
legend("$\vert \tilde{\alpha_N}^{\Delta} - \alpha \vert$", ...
       "interpreter", "latex", "Location", "Best", "FontSize", 10);
title("Convergence of the estimator $\tilde{\alpha_N}^{\Delta}$",  "interpreter", "latex");
saveas(gcf, "plot_Q12", "epsc");