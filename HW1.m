clear all; close all; clc;

%% GENERATE THE PROCESS x(k), 1 REALIZATION

Nsamples=800;
% Frequencies of the exponentials
f1=0.17;
f2=0.78;
% Generate the white noise (2 components)
sigmaw=0.1;
% Real part
wi=sigmaw*randn(Nsamples,1);
% Imaginary part
wq=sigmaw*randn(Nsamples,1);
% Generate the initial phases
phi1=2*pi*rand(1);
phi2=2*pi*rand(1);

xi=zeros(Nsamples,1);
xq=zeros(Nsamples,1);
for k=1:Nsamples
    xi(k)=cos(2*pi*f1*k+phi1)+cos(2*pi*f2*k+phi2)+wi(k);
    xq(k)=sin(2*pi*f1*k+phi1)+sin(2*pi*f2*k+phi2)+wq(k);
end

% Complex r.p. x(k), 800 samples
x=xi+j*xq;

%% SPECTRAL ANALYSIS

% Autocorrelation (unbiased estimate)
[rx]=autocorrelation_Unb(x);
L=floor(Nsamples/3);   %L should be lower than the length of the r.p. because of the high variance when n approaches K
rx=rx(1:L);

% Blackman-Tukey correlogram 
% The length of 2*L+1 is because of page 86 note 24
w_rect=window(@rectwin,2*L+1);
w_hamming=window(@hamming,2*L+1);    
Pbt1=correlogram(x, w_hamming, rx, L);
Pbt2=correlogram(x, w_rect, rx, L);

figure();
subplot(2,1,1);
plot(1/Nsamples:1/Nsamples:1,10*log10(abs(Pbt1)));
title('Correlogram - Hamming');
ylabel('Amplitude (dB)');
xlabel('f');
subplot(2,1,2);
plot(1/Nsamples:1/Nsamples:1,10*log10(abs(Pbt2)));
title('Correlogram - rect');
xlabel('f');
ylabel('Amplitude (dB)');
ylim([-15 30]);

% Periodogram
X=fft(x);
Pper=(1/Nsamples)*(abs(X)).^2;
figure,
plot(1/Nsamples:1/Nsamples:1,10*log10(Pper))
title('Periodogram estimate of the PSD')
xlabel('f')
ylabel('Amplitude (dB)')
ylim([-15 30])

% Welch periodogram
w_welch=window(@hamming,100);
Welch_P = welchPSD(x, w_welch, 20);
figure('Name','Welch periodogram');
plot(1/Nsamples:1/Nsamples:1,10*log10(Welch_P))
title('Welch periodogram estimate of the PSD')
xlabel('f')
ylabel('Amplitude (dB)')
ylim([-15 30])

% Analytical PSD: compute the transform of rx(n) on paper and plot it
% according to the requirements

%% Optimal choice of N

err_mean = zeros(1, 10);

for i = 1:20
    % Find c_opt using N coefficients
    c = findAR(i, rx);
    % Error is taken on 100 different samples
    e = zeros(1,100);
    for k = 1:100
        if (k < i + 1)
            % Input vector of length N
            x_in = flipud([zeros(i-k+1,1); x(1:k-1)]);
            % For k = 1 x(1:0) is an empty matrix
            y_k = x_in.'*c(1:i);
        else
            % Revert vector to obtain values from k-1 to k-N
            x_in = flipud(x((k-i):(k-1)));
            y_k = x_in.'*c(1:i);
        end
        % Computing the error, d(k) = x(k)
        e_k = x(k) - y_k;
        e(k) = e_k;
    end
    err_mean(i) = abs(sum(e))/length(e); 
end

figure('Name', 'Error over N');
title('Error in function of N');
plot(1:20, abs(err_mean));
xlim([1 20]);
xlabel('N'); ylabel('Error');

% Setting N as argmin of error
[minerr, N]  = min(err_mean);
c_opt = -findAR(N, rx);

%% AR model
% Coefficients of Wiener filter
[a, s_white, d]=findAR(N, rx);
[H_w, omega] = freqz(1, [1; a], Nsamples, 'whole');
figure('Name','AR model estimate of the PSD');
title('AR model estimate of the PSD');
plot(omega/(2*pi), 10*log10(s_white*abs(H_w).^2));
xlabel('f');
ylabel('Amplitude (dB)');
ylim([-15 40]);

%% Final spectral plot
figure('Name', 'Spectral Analysis');
title('Spectral analysis');
hold on;
plot((1:Nsamples)/Nsamples, 10*log10(Welch_P), 'r-.')
plot((1:Nsamples)/Nsamples, 10*log10(abs(Pbt1)), 'Color', 'b')
plot((1:Nsamples)/Nsamples, 10*log10(Pper), 'g:')
plot(omega/(2*pi), 10*log10(s_white*abs(H_w).^2), 'Color', 'm');
legend('Welch', 'Correlogram', 'Periodogram', ['AR(' int2str(N) ')'], 'Location', 'SouthWest');
hold off;

[H, www] = freqz([1; a], 1, Nsamples, 'whole');
figure('Name', 'Z-plane for error predictor A(z)');
title('Z-plane for error predictor A(z)');
zplane([1;a]);

%% Least Mean Squares Filter

% Length
K = L;
% Max number of iterations
max_iter = 800;
% Coefficients and error initialization
c = zeros(N, max_iter + 1);
e = zeros(1, max_iter);
% Choice of parameter mu
mu_tilde = 1;
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

for index = 1:N
    figure('Name', ['Coefficient of index ' int2str(index)]);
    subplot(2, 1, 1)
    plot(1:max_iter+1, real(c(index, :)));
    hold on;
    plot([1, max_iter+1], -real(a(index))*[1 1]);
    title(['Real part of c' int2str(index) ' and c_{opt}' int2str(index)]);
    legend(['c' int2str(index)], ['a' int2str(index)]);
    xlim([0 max_iter]);
    
    subplot(2, 1, 2);
    plot(1:max_iter+1, imag(c(index, :)));
    hold on;
    plot([1,max_iter+1], -imag(a(index))*[1 1]);
    title(['Imaginary part of c' int2str(index) ' and c_{opt}' int2str(index)]);
    legend(['c' int2str(index)], ['a' int2str(index)]);
    xlim([0 max_iter]);
end

figure('Name', 'Error function');
title('Error function at each iteration');
plot(1:max_iter, 10*log10(abs(e).^2), 1:max_iter, 10*log10(s_white)*ones(1, max_iter));

%% 3 PREDICTOR
[copt, Jmin]=predictor(rx, N);
t=20;
Jvect=zeros(t,1);

for i=1:length(Jvect)
    [c_it, J_it]=predictor(rx, i);
    Jvect(i)=J_it;
end

figure, plot([1:t],Jvect)
coeff=[1; copt];
A = tf([1 copt.'], 1,1);
figure, pzmap(A)
