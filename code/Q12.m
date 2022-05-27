%% parameters
clear 
close all
clc

rng(0)
T = 1e3; % final time
M = 1e4; % averages
h = 2^(-10); % step size for EM
Deltas = 2.^(-(0:7)); 
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition

%% estimate of \alpha using Q11)
estimates_alpha = zeros(length(Deltas), 1);
W = simulate_brownian(0, T, T/h-1);
X = euler_maruyama(x0, alpha, sigma, h, W);
for i=1:length(Deltas)
    X_observed = X(1:Deltas(i)/h:end);
    estimates_alpha(i) = estimator_alpha_Q11(Deltas(i), X_observed);
end

%% plot of estimated \alpha
figure()
semilogx(1./Deltas', estimates_alpha, 'Linewidth', 1.5);
hold on
semilogx(1./Deltas', ones(length(Deltas), 1), '--r', 'DisplayName', "$\alpha = " + num2str(alpha) + "$");
legend("$\widetilde{\alpha}_N^{\Delta}$", "$\alpha = $" + num2str(alpha), ...
       "interpreter", "latex", "Location", "Best", "FontSize", 20);
xlabel("$1/\Delta$", "interpreter", "latex");
ylim([min(estimates_alpha)-0.1 max([sigma+0.05, max(estimates_alpha)])]);
title("Behaviour of $\widetilde{\alpha}_N^{\Delta}$ against $\Delta$",  "interpreter", "latex", "FontSize", 13);
saveas(gcf, "plot_Q12", "epsc");