clc; close all; clear global; clearvars;

%% Least Mean Squares Estimation

% Load one realization
load('realizations.mat','x');

% Max number of iterations
max_iter = 800;
% Set parameters
L = floor(size(x,1)/3);
N = 3;
K = L;

% Error vector initialization
e = zeros(size(x,2), max_iter);

for i=1:size(x,2)
    % Autocorrelation
    rx = autocorrelation_Unb(x(:,i));
    rx = rx(1:L);
    [a, s_white] = findAR(N, rx);
    % Coefficients initialization
    c = zeros(N, max_iter + 1);
    % Choice of parameter mu
    mu_tilde = 1;
    mu = mu_tilde/(rx(1)*N);

    % Center signal around its mean
    z = x(:,i) - mean(x(:,i));

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
                e(i,k) = e_k;
                % Update step
                c(:, k+1) = c(:, k) + mu*e_k*conj(x_in);
        end
end