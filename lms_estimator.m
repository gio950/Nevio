clc; close all; clear global; clearvars;

%% Least Mean Squares Estimation

% Load one realization
load('inputsignal01.mat','x');

% Set parameters
L = floor(length(x)/3);
N = 2;
rx = autocorrelation_Unb(x);
rx = rx(1:L);
[a, s_white] = findAR(N, rx);
K = L;

% Max number of iterations
max_iter = 800;

% Coefficients and error initialization
c = zeros(N, max_iter + 1);
e = zeros(1, max_iter);

% Choice of parameter mu
mu_tilde = 0.05;
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

for index = 1:N
    figure('Name', ['Coefficient of index ' int2str(index)]);
    subplot(2, 1, 1)
    plot(1:max_iter+1, real(c(index, :)));
    hold on;
    plot([1, max_iter+1], -real(a(index))*[1 1]);
    hold on;
    plot(1:max_iter+1, real(c_mean(index, :)));
    title(['Real part of c' int2str(index) ' and c_{opt}' int2str(index)]);
    legend(['c' int2str(index)], ['c_{est}' int2str(index)], ['average' int2str(index)]);
    xlim([0 max_iter]);
    
    
    subplot(2, 1, 2);
    plot(1:max_iter+1, imag(c(index, :)));
    hold on;
    plot([1,max_iter+1], -imag(a(index))*[1 1]);
    hold on;
    plot(1:max_iter+1, imag(c_mean(index, :)));
    title(['Imaginary part of c' int2str(index) ' and c_{opt}' int2str(index)]);
    legend(['c' int2str(index)], ['c_{est}' int2str(index)],['average' int2str(index)]);
    xlim([0 max_iter]);
    xlabel('Number of iterations')
end

figure('Name', 'Error function');
plot(1:max_iter, 10*log10(abs(e).^2), 1:max_iter, 10*log10(s_white)*ones(1, max_iter), 1:max_iter, 10*log10(mean_error));

title('Error function at each iteration');
legend('|e(k)|^2', '\sigma_w^2','J_{min}');
ylim([-15 10])
xlabel('k')
ylabel('|e(k)|^2, J_{min}, \sigma_w^2')