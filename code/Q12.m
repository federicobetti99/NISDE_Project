%% parameters
rng(0)
T = 1e3; % final time
M = 1e4; % averages
h = 2^(-10); % step size for EM
Deltas = 2.^(-(0:7)); 
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition

%% estimate of \alpha using Q11)
error_estimates_alpha = zeros(length(Deltas), 1);
W = simulate_brownian(0, T, T/h-1);
X = euler_maruyama(x0, alpha, sigma, h, W);
for i=1:length(Deltas)
    X_observed = X(1:Deltas(i)/h:end);
    error_estimates_alpha(i) = abs(estimator_alpha_Q11(Deltas(i), X_observed)-alpha);
end

%% plot of estimated \alpha
figure()
semilogx(Deltas, error_estimates_alpha, 'Linewidth', 1.0);
legend("$\vert \tilde{\alpha_N}^{\Delta} - \alpha \vert$", ...
       "interpreter", "latex", "Location", "Best", "FontSize", 10);
xlabel("$\Delta$", "interpreter", "latex");
title("Convergence of the estimator $\tilde{\alpha_N}^{\Delta}$",  "interpreter", "latex", "FontSize", 15);
saveas(gcf, "plot_Q12", "epsc");