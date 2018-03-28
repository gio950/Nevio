close all; clear global; clearvars; clc;

%% GENERATE THE PROCESS x(k), 1 REALIZATION

Nsamples=800;
% Frequencies of the exponentials
f1=0.17;
f2=0.78;
% Generate the white noise (2 components)
sigmaw=2;
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
x=xi+1i*xq;

%
[rx]=autocorrelation_Unb(x);
L=floor(Nsamples/5);%L should be lower than the length of the r.p. because of the high variance when n approaches K
rx=rx(1:L);
N = 10;
[copt, Jmin ]=predictor(rx, N);
coeff=[1; copt];
[a, s_white, d]=findAR(N, rx);
zplane([1; a].', 1);
[H_w, omega] = freqz(1, [1; a], Nsamples, 'whole');
figure(), plot(omega/(2*pi), 10*log10(s_white*(abs(H_w)).^2), 'Color', 'm');

% save('aaaaa.mat', 'x');