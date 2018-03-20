clearvars
close all
clc

%% GENERATE THE PROCESS x(k), 1 REALIZATION

Nsamples=800;
%frequencies of the exponentials
f1=0.17;
f2=0.78;
%generate the white noise (2 components)
sigmaw=2;
%real part
wi=sigmaw*randn(Nsamples,1);
%imaginary part
wq=sigmaw*randn(Nsamples,1);
%generate the initial phases
phi1=2*pi*rand(1);
phi2=2*pi*rand(1);

xi=zeros(Nsamples,1);
xq=zeros(Nsamples,1);
for k=1:Nsamples
    xi(k)=cos(2*pi*f1*k+phi1)+cos(2*pi*f2*k+phi2)+wi(k);
    xq(k)=sin(2*pi*f1*k+phi1)+sin(2*pi*f2*k+phi2)+wq(k);
end

%complex r.p. x(k) 800 samples
x=xi+j*xq;

%% SPECTRAL ANALYSIS

%autocorrelation (unbiased estimate)
[rx]=autocorrelation_Unb(x);
L=floor(Nsamples/3);   %L should be lower than the length of the r.p. because of the high variance when n approaches K
rx=rx(1:L);

%Blackman-Tukey correlogram 
%the length of 2*L+1 is because of page 86 note 24
w_rect=window(@rectwin,2*L+1);
w_hamming=window(@hamming,2*L+1);    
Pbt1=correlogram(x, w_hamming, rx, L);
Pbt2=correlogram(x, w_rect, rx, L);

figure,
subplot(211)
plot(1/Nsamples:1/Nsamples:1,10*log10(abs(Pbt1)))
title('Correlogram - Hamming')
ylabel('Amplitude (dB)')
xlabel('f')
subplot(212)
plot(1/Nsamples:1/Nsamples:1,10*log10(abs(Pbt2)))
title('Correlogram - rect')
xlabel('f')
ylabel('Amplitude (dB)')
ylim([-15 30])

%Periodogram

X=fft(x);
Pper=(1/Nsamples)*(abs(X)).^2;
figure,
plot(1/Nsamples:1/Nsamples:1,10*log10(Pper))
title('Periodogram estimate of the PSD')
xlabel('f')
ylabel('Amplitude (dB)')
ylim([-15 30])

% Welch periodogram
Welch_P = welchPSD(x, w_hamming, 25);
figure('Name','Welch periodogram');
plot(1/Nsamples:1/Nsamples:1,10*log10(Welch_P))
title('Welch periodogram estimate of the PSD')
xlabel('f')
ylabel('Amplitude (dB)')
ylim([-15 30])

%analytical PSD: compute the transform of rx(n) on paper and plot it
%according to the requirements
 
 
%% AR model
% Coefficients of Wiener filter
N=40;
[a, s_white, d]=findAR(N, rx);
[H, omega] = freqz(1, [1; a], Nsamples, 'whole');
figure,
plot(omega/(2*pi), 10*log10(s_white*abs(H).^2))
title('AR model estimate of the PSD')
xlabel('f')
ylabel('Amplitude (dB)')
ylim([-15 40])

%% final plot
figure, hold on
plot((1:Nsamples)/Nsamples, 10*log10(Welch_P), 'r-.')
plot((1:Nsamples)/Nsamples, 10*log10(abs(Pbt1)), 'Color', 'b')
plot((1:Nsamples)/Nsamples, 10*log10(Pper), 'c:')
plot(omega/(2*pi), 10*log10(s_white*abs(H).^2), 'Color', 'm');
legend('Welch', 'Correlogram', 'Periodogram', ['AR(' int2str(N) ')'], 'Location', 'SouthWest')
hold off
title('Spectral analysis')


[H, www] = freqz([1; a], 1, Nsamples, 'whole');
figure('Name', 'Z-plane for error predictor A(z)');
title('Z-plane for error predictor A(z)');
zplane([1;a]);

%% Least Mean Squares Filter

% Length
K = L;
% Look up in the sky and set N
N = 4;
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
        z_k_1 = flipud([zeros(N - k + 1, 1); z(1:k - 1)]);
        % For k = 1 z(1:0) is an empty matrix
        y_k = z_k_1.'*c(:, k);
    else
        % Revert vector to obtain values from k-1 to k-N
        z_k_1 = flipud(z((k - N):(k-1)));
        y_k = z_k_1.'*c(:, k);
    end
    % d(k) is input signal z(k)
    e_k = z(k) - y_k;
    e(k) = e_k;
    % Update step
    c(:, k+1) = c(:, k) + mu*e_k*conj(z_k_1);
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