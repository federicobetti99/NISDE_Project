%% parameters
clear
close all
clc

rng('default')
T = 1e3; % final time
M = 1e4; % averages
h = 2^(-10); % step size for EM
Deltas = 2.^(-(0:7)); 
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition

%% compute estimates of \alpha and \sigma
estimates_alpha = zeros(length(Deltas), 1);
estimates_sigma = zeros(length(Deltas), 1);
X = euler_maruyama(x0, alpha, sigma, h, T); %% solve with EM
for i=1:length(Deltas)
    X_observed = X(1:(Deltas(i)/h):end); %% get only the observed values
    estimates_alpha(i) = estimator_alpha_Q7(Deltas(i), X_observed); %% compute estimate for \alpha
    estimates_sigma(i) = estimator_sigma_Q7(Deltas(i), X_observed); %% compute estimate for \sigma
end

%% plot of estimated \alpha and \sigma against the true ones
figure()
semilogx(Deltas', estimates_alpha, '-o', 'Linewidth', 1.5);
hold on
semilogx(Deltas', ones(length(Deltas), 1), '--r', 'DisplayName', "$\alpha = " + num2str(alpha) + "$");
xlabel("$\Delta$", "interpreter", "latex");
ylim([min(estimates_alpha)-0.1 max([sigma+0.05, max(estimates_sigma)])]);
legend("$\widehat{\alpha}_N^{\Delta}$", "$\alpha = $" + num2str(alpha), "location", "best", "FontSize", 20, "interpreter", "latex");
title("Behaviour of $\widehat{\alpha}_N^{\Delta}$ against $\Delta$", "interpreter", "latex", "FontSize", 13);
saveas(gcf, "plot_Q7_alpha", "epsc");

figure()
semilogx(Deltas', estimates_sigma, '-o', 'Linewidth', 1.5);
hold on
semilogx(Deltas', ones(length(Deltas), 1), '--r', 'DisplayName', "$\sigma = " + num2str(sigma) +"$");
ylim([min(estimates_sigma)-0.1 max([sigma+0.05, max(estimates_sigma)])]);
xlabel("$\Delta$", "interpreter", "latex");
legend("$\widehat{\sigma}_N^{\Delta}$", "$\sigma = $" + num2str(sigma), "location", "best", "FontSize", 20, "interpreter", "latex");
title("Behaviour of $\widehat{\sigma}_N^{\Delta}$ against $\Delta$",  "interpreter", "latex", "FontSize", 13);
saveas(gcf, "plot_Q7_sigma", "epsc");