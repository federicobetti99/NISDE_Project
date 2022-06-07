%% parameters
clear
close all
clc

rng('default')
T = 1e3; % final time
M = 1e4; % averages
h = 0.01; % step size for EM
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition
pho_infty = @(x) 1/sqrt(2*pi*sigma/alpha) * exp(-x.^2 / (2*sigma/alpha)); % invariant density

%% calculate empirical averages
timegrid = linspace(0, T, T/h);
mean_endpoint = zeros(M, 1);
for k = 1:M
    X = euler_maruyama(x0, alpha, sigma, h, T);
    mean_endpoint(k) = X(end);
end

%% plot histogram
x = linspace(-5, 5, 10000);
figure()
histogram(mean_endpoint, 'normalization', 'pdf');
hold on
plot(x, pho_infty(x), 'LineWidth', 1.5);
ylim([0 1/sqrt(2*pi*sigma/alpha)+0.1]);
legend('', "\rho_{\infty}(x)", 'Fontsize', 15);
title("Histogram of $\{X^{(m)}(T)\}_{m = 1}^{M}$, where $T = 10^3$ and $M = 10^4$", "interpreter", "latex", "FontSize", 13);
saveas(gcf, "plot_Q3", "epsc");
