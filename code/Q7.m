%% parameters
rng(0)
T = 1e3; % final time
M = 1e4; % averages
h = 2^(-10); % step size for EM
Deltas = 2.^(-(0:7)); 
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition

%% compute estimators of \alpha and \sigma
error_estimates_alpha = zeros(length(Deltas), 1);
error_estimates_sigma = zeros(length(Deltas), 1);
W = simulate_brownian(0, T, T/h-1); %% realization of a BM
X = euler_maruyama(x0, alpha, sigma, h, W); %% solve with EM
for i=1:length(Deltas)
    X_observed = X(1:(Deltas(i)/h):end); %% get only the observed values
    error_estimates_alpha(i) = abs(estimator_alpha(Deltas(i), X_observed)-alpha); %% compute estimate for \alpha
    error_estimates_sigma(i) = abs(estimator_sigma(Deltas(i), X_observed)-sigma); %% compute estimate for \sigma
end

%% plot of estimated \alpha and \sigma against the true ones
figure()
semilogx(Deltas, error_estimates_alpha, 'Linewidth', 1.0);
hold on
semilogx(Deltas, error_estimates_sigma, 'Linewidth', 1.0);
xlabel("$\Delta$", "interpreter", "latex");
legend("$\vert \hat{\alpha_N}^{\Delta} - \alpha \vert$", "$ \vert \hat{\sigma_N}^{\Delta} - \sigma \vert$", ...
       "interpreter", "latex", "Location", "Best", "FontSize", 10);
title("Convergence of the estimators $\hat{\alpha_N}^{\Delta}$ and $\hat{\sigma_N}^{\Delta}$",  "interpreter", "latex", "FontSize", 15);
saveas(gcf, "plot_Q7", "epsc");

