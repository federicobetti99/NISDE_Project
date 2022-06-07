%% parameters
clear
close all
clc

rng('default')
T = 5 * 1e3; % final time
h = 0.01; % step size for EM
Delta = 1; % sampling rate
alpha = 1; % drift coefficient
sigma = 1; % diffusion coefficient
x0 = 2; % initial condition

%% simulation
% compute numerical solution
X = euler_maruyama(x0, alpha, sigma, h, T);

% Sampling the solution with sampling rate delta
N = T/Delta;
skip = (length(X) - 1) / N;
X = X(1:skip:end);

%% solve non-linear system
Ns = 2:N;   % Numbers of available observations
x0 = [0.5; 0.5];
alpha_est = zeros(length(Ns), 1);
sigma_est = zeros(length(Ns), 1);

for i = 1:length(Ns)
    n = Ns(i);
    X_n = X(1:n+1);
    fun = @(param) G_2D_Q15(param, X_n, Delta);
    options = optimset('Display','off');
    x = fsolve(fun,x0,options); % Compute estimates
    alpha_est(i) = x(1);
    sigma_est(i) = x(2);
end


%% plot of estimated alpha_n and sigma_n against the number of available observations n
fig = figure();
plot(linspace(2, N, N-1), sigma_est, 'DisplayName', '$\widetilde{\sigma}_n^{\Delta}$', LineWidth=1)
hold on
plot(linspace(2, N, N-1), ones(N-1, 1), '--r', 'DisplayName', "$\sigma = " + num2str(sigma) +"$")
ylim([-inf max([max(sigma_est), 1]) + 0.02])
xlabel('$n$', 'Interpreter','latex', 'FontSize',16)
legend('Interpreter','latex', 'Location', 'best', 'FontSize',16)
saveas(fig, 'Q15_sigma', 'epsc')

fig = figure();
plot(linspace(2, N, N-1), alpha_est, 'DisplayName', '$\widetilde{\alpha}_n^{\Delta}$', LineWidth=1)
hold on
plot(linspace(2, N, N-1), ones(N-1,1), '--r', 'DisplayName', "$\alpha = " + num2str(alpha) +"$")
ylim([-inf max([max(alpha_est), 1]) + 0.02])
xlabel('$n$', 'Interpreter','latex', 'FontSize',16)
legend('Interpreter', 'latex', 'Location','best', 'FontSize',16)
saveas(fig, 'Q15_alpha', 'epsc')