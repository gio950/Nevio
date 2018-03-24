clc; close all; clear global; clearvars;

%% Least Mean Squares Estimation

% Load one realization
load('realizations.mat','x');

% Max number of iterations
max_iter = 800;
% Set parameters
L = floor(size(x,1)/3);
N = 2;
K = L;

% Error vector initialization
e = zeros(size(x,2), max_iter);
C1=zeros(300,801);
C2=zeros(300,801);
for i=1:size(x,2)
    % Autocorrelation
    rx = autocorrelation_Unb(x(:,i));
    rx = rx(1:L);
    [a, s_white] = findAR(N, rx);
    % Coefficients initialization
    c = zeros(N, max_iter + 1);
    % Choice of parameter mu
    mu_tilde = 0.1;
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
        c1=c(1,:);
        c2=c(2,:);
        C1(i,:)=c1;
        C2(i,:)=c2;
end

c_mean(1,:)=mean(C1);
c_mean(2,:)=mean(C2);
% Mean error is computed over 300 errors for the same k

%mean_error = zeros(1, max_iter);
mean_error = mean(abs(e.^2));
%for k=1:max_iter
    %mean_error(k) = sum(e(:,k))/size(e,1);
%end


for index = 1:N
    figure('Name', ['Coefficient of index ' int2str(index)]);
    subplot(2, 1, 1)
    plot(1:max_iter+1, real(c_mean(index, :)));
    hold on;
    plot([1, max_iter+1], -real(a(index))*[1 1]);
    title(['Real part of c_{mean}' int2str(index) ' and c_{opt}' int2str(index)]);
    legend(['c' int2str(index)], ['a' int2str(index)]);
    xlim([0 max_iter]);
    
    subplot(2, 1, 2);
    plot(1:max_iter+1, imag(c_mean(index, :)));
    hold on;
    plot([1,max_iter+1], -imag(a(index))*[1 1]);
    title(['Imaginary part of c_{mean}' int2str(index) ' and c_{opt}' int2str(index)]);
    legend(['c' int2str(index)], ['a' int2str(index)]);
    xlim([0 max_iter]);
end

figure('Name','Mean squared error');
plot(1:max_iter,10*log10(mean_error), 1:max_iter, 10*log10(s_white)*ones(1, max_iter));
title('Mean squared error over iterations');
xlabel('k (iterations)'); ylabel('Mean Squared Error (dB)');


save('Jmin.mat', 'mean_error');
save('avg_coeff.mat', 'c_mean');