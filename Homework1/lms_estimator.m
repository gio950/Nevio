clc; close all; clear global; clearvars;

%% Least Mean Squares Estimation

% Load one realization
load('inputsignal01.mat','x');

% Set parameters
L = floor(length(x)/5);
N = 2;
rx = autocorrelation_Unb(x);
rx = rx(1:L);
[a, s_white] = findAR(N, rx);
K = L;

% Autocorrelation Matrix
col = rx(1:N);
row = conj(col);
R=toeplitz(col, row);
mu_opt = 2/sum(eig(R));

% Max number of iterations
max_iter = 800;

% Coefficients and error initialization
c = zeros(N, max_iter + 1);
e = zeros(1, max_iter);

% Choice of parameter mu
mu_tilde = 0.06;
mu = mu_tilde/(rx(1)*N);

% Center signal around its mean
z = x - mean(x);

for k = 1:max_iter
    if (k < N + 1)
        % Input vector of length N
        x_in = flipud([zeros(N - k + 1, 1); z(1:k - 1)]);
        % For k = 1 z(1:0) is an empty matrix
        y_k = x_in.'*c(:, k);
    else
        % Revert vector to obtain values from k-1 to k-N
        x_in = flipud(z((k - N):(k-1)));
        y_k = x_in.'*c(:, k);
    end
    % d(k) is input signal z(k)
    e_k = z(k) - y_k;
    e(k) = e_k;
    % Update step
    c(:, k+1) = c(:, k) + mu*e_k*conj(x_in);
end

load('Jmin.mat', 'mean_error');
load('avg_coeff.mat', 'c_mean');

% set(0,'defaultTextInterpreter','latex')          % to use LaTeX format
% set(gca,'FontSize',10);
%% Plot of Real Part
figure()
line(1:max_iter+1,real(c(1, :)),'Color','r');
line(1:max_iter+1,-real(a(1))*ones([1,max_iter+1]),'Color','k','LineStyle','--');
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k'
ax1_pos = ax1.Position; % position of first axes
xlabel('Number of iterations')
ylim([0 1]);
xlim([1 800]);
ylabel('$c_1$');
legend('c_{1,lms}');

ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
line(1:max_iter+1,real(c(2, :)),'Color','b');
line(1:max_iter+1,-real(a(2))*ones([1,max_iter+1]),'Color','k','LineStyle','--');
ylim([-1 0]);
xlim([1 800]);
ylabel('$c_2$');
legend('c_{2,lms}');

%% Plot of Imaginary Part
figure()
line(1:max_iter+1,imag(c(1, :)),'Color','r');
line(1:max_iter+1,-imag(a(1))*ones([1,max_iter+1]),'Color','k','LineStyle','--');
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k'
ax1_pos = ax1.Position; % position of first axes
xlabel('Number of iterations')
ylim([-0.15 0.1]);
xlim([1 800]);
ylabel('$c_1$');
legend('c_{1,lms}');

ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
line(1:max_iter+1,imag(c(2, :)),'Color','b');
line(1:max_iter+1,-imag(a(2))*ones([1,max_iter+1]),'Color','k','LineStyle','--');
ylim([-0.05 0.39]);
xlim([1 800]);
ylabel('$c_2$');
legend('c_{2,lms}');

%% Plot of Jmin
figure('Name', 'Error function');
plot(1:max_iter, 10*log10(abs(e).^2)); hold on
plot(1:max_iter, 10*log10(mean_error)); hold on
plot(1:max_iter, 10*log10(s_white)*ones(1, max_iter)','r--','LineWidth',2);
grid on
title('Error function at each iteration','FontSize',15);
legend('|e(k)|^2','J(k)','J_{min}');
ylim([-15 10])
xlabel('k')
ylabel('$|e(k)|^2$,$J(k)$, $J_{min}$')


set(0,'defaultTextInterpreter','latex')          % to use LaTeX format
set(gca,'FontSize',20);

