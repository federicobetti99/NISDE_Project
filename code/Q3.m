%% parameters
T = 1e3; % final time
M = 1e4; % averages
h = 0.01; % step size for EM
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition
var = @(t) sigma/alpha * (1-exp(-2*alpha*t)); %% variance of X(t)
pho_infty = @(x) 1/sqrt(2*pi*sigma/alpha) * exp(-x.^2 / (2*sigma/alpha)); % invariant density

%% calculate empirical averages
timegrid = linspace(0, T, T/h);
mean_endpoint = zeros(M, 1);
for k = 1:M
    W = simulate_brownian(0, T, length(timegrid)-1);
    X = euler_maruyama(x0, alpha, sigma, h, W);
    mean_endpoint(k) = X(end);
end

%% example of solution
figure()
plot(timegrid(1:100), X(1:100))
title("Sample path of the Ornstein–Uhlenbeck process", "FontSize", 10);
saveas(gcf, "example_sol", "epsc");


%% plot histograms
x = linspace(-5, 5, 10000);
figure()
histogram(mean_endpoint, 'normalization', 'pdf');
hold on
plot(x, pho_infty(x), 'r', 'LineWidth', 1.5);
ylim([0 1/sqrt(2*pi*sigma/alpha)+0.1]);
legend('', "\rho_{\infty}(x)", 'Fontsize', 15);
title("Empirical and theoretical distribution of X(T), $T = 10^3$", "interpreter", "latex", "FontSize", 13);
saveas(gcf, "plot_Q3", "epsc");
