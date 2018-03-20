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
