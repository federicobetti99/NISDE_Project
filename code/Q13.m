%% parameters
clear 
close all
clc

rng(0)
T = 1e3; % final time
M = 1e4; % averages
h = 0.01; % step size for EM
Delta = 1; % sampling rate
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition
Sigma = 1/Delta^2 * (exp(2*alpha*Delta)-1);
pho_sigma = @(x) 1/sqrt(2*pi*Sigma) * exp(-x.^2/(2*Sigma));

%% measuring estimates
clt_values = zeros(M, 1);
spacegrid = linspace(-10, 10);
N = T/Delta;
for k=1:M
    W = simulate_brownian(0, T, T/h-1);
    X = euler_maruyama(x0, alpha, sigma, h, W);
    X_observed = X(1:Delta/h:end);
    clt_values(k) = sqrt(N) * (estimator_alpha_Q11(Delta, X_observed)-alpha);
end

%% plot
figure()
histogram(clt_values, 'normalization', 'pdf');
hold on
plot(spacegrid, pho_sigma(spacegrid), 'LineWidth', 2);
legend("", "$\mathcal{N}(0, \Sigma)$", "interpreter", "latex", "Fontsize", 20);
xlim([-10 10]);
ylim([0 1/sqrt(2*pi*Sigma)+0.05]);
saveas(gcf, "plot_Q13", "epsc");