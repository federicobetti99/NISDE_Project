%% parameters
clear 
close all
clc

rng('default')
T = 1e3; % final time
M = 1e4; % averages
h = 0.01; % step size for EM
Delta = 1; % sampling rate
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition
Sigma = 1/Delta^2 * (exp(2*alpha*Delta)-1); % variance of the estimator
pho_sigma = @(x) 1/sqrt(2*pi*Sigma) * exp(-x.^2/(2*Sigma)); % asymptotic density

%% estimating \alpha
clt_values = zeros(M, 1);
spacegrid = linspace(-10, 10);
N = T/Delta;
for k=1:M
    X = euler_maruyama(x0, alpha, sigma, h, T);
    X_observed = X(1:Delta/h:end);
    clt_values(k) = sqrt(N) * (estimator_alpha_Q12(Delta, X_observed)-alpha);
end

%% plot of the distribution
figure()
histogram(clt_values, 'normalization', 'pdf');
hold on
plot(spacegrid, pho_sigma(spacegrid), 'LineWidth', 2);
legend("", "$\rho_{\mu_{\alpha, \Delta}}(x)$", "interpreter", "latex", "Fontsize", 20);
xlim([-10 10]);
ylim([0 1/sqrt(2*pi*Sigma)+0.05]);
title("Histogram of $\sqrt{N} \left(\{\tilde{\alpha}_{N}^{\Delta, (m)}\}_{m=1}^{M} - \alpha \right)$, $M = 10^4$",  "FontSize", 13, "interpreter", "latex")
saveas(gcf, "plot_Q13", "epsc");