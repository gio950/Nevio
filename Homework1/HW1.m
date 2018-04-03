
clc; 
close all; 
clear global; clearvars;

%% LOAD 1 REALIZATION OF THE PROCESS

load('inputsignal01.mat', 'x');
Nsamples=length(x);

%% SPECTRAL ANALYSIS

% Autocorrelation (unbiased estimate)
[rx]=autocorrelation_Unb(x);
L=floor(Nsamples/5);% L should be lower than the length of the r.p. because of the high variance when n approaches K
rx=rx(1:L);

% Blackman-Tukey correlogram 
% The length of 2*L+1 is because of page 86 note 24
w_rect=window(@rectwin,2*L+1);
w_hamming=window(@hamming,2*L+1); 
w_bartlett=window(@bartlett,2*L+1);
Pbt3=correlogram(x, w_hamming, rx, L);
Pbt2=correlogram(x, w_rect, rx, L);
Pbt1=correlogram(x, w_bartlett, rx, L);    % use Pbt3 for sigma=2

% comparison of different window types
% figure();
% hold on
% plot(1/Nsamples:1/Nsamples:1,10*log10(abs(Pbt1)),'g'); 
% plot(1/Nsamples:1/Nsamples:1,10*log10(abs(Pbt2)));
% plot(1/Nsamples:1/Nsamples:1,10*log10(abs(Pbt3)),'r','Linewidth',2);
% title('Correlogram vs window type');
% xlabel('f');
% ylabel('Amplitude (dB)');
% ylim([-15 30]);
% legend('Hamming','Rect','Bartlett');


% Periodogram
X=fft(x);
Pper=(1/Nsamples)*(abs(X)).^2;

%{
figure,
plot(1/Nsamples:1/Nsamples:1,10*log10(Pper))
title('Periodogram estimate of the PSD')
xlabel('f')
ylabel('Amplitude (dB)')
ylim([-15 30])
%}

% Welch periodogram
S=70;   %overlap
D=100;   %window length
w_welch=window(@hamming,D);
[Welch_P, Ns] = welchPSD(x, w_welch, S);
var_Welch=Welch_P.^2/Ns;

% figure('Name','Welch periodogram');
% plot(1/Nsamples:1/Nsamples:1,10*log10(Welch_P))
% title('Welch periodogram estimate of the PSD')
% xlabel('f')
% ylabel('Amplitude (dB)')
% ylim([-15 30])

% Comparison of different S,D
% S = [20 50 70 200];
% D = [50 100 150 400];
% for i=1:4
%     w_welch=window(@hamming,D(i));
%     [Welch_P(:,i), Ns] = welchPSD(x, w_welch, S(i));
%     var_Welch=Welch_P(:,i).^2/Ns;
% end
% figure('Name','Welch periodogram as a function of D and S');
% plot(1/Nsamples:1/Nsamples:1,10*log10(Welch_P))
% title('Welch periodogram estimate of the PSD')
% xlabel('f')
% ylabel('Amplitude (dB)')
% ylim([0 25])
% grid
% legend(['D = ' int2str(D(1)) ' and S = ' int2str(S(1)) ], ['D = ' int2str(D(2)) ' and S = ' int2str(S(2)) ],['D = ' int2str(D(3)) ' and S = ' int2str(S(3)) ],['D = ' int2str(D(4)) ' and S = ' int2str(S(4)) ]);
%     


% Analytical PSD: compute the transform of rx(n) on paper and plot it
% according to the requirements
%sigmaw=0.1
b = zeros(1,800);
for i=1:length(b)
    b(i) = 10*log10(0.1);
end

% 30 is random choice just to see the plot
b(ceil(0.17*800)) = 10*log10(Nsamples);
b(ceil(0.78*800)) = 0.8*10*log10(Nsamples);

%sigmaw=2

% b = zeros(1,800);
% for i=1:length(b)
%     b(i) = 10*log10(2);
% end
% 
% % 30 is random choice just to see the plot
% b(ceil(0.17*800)) = 10*log10(Nsamples);
% b(ceil(0.78*800)) = 0.8*10*log10(Nsamples);

%% Choice of N
N =5;
[copt, Jmin, det_R ]=predictor(rx, N);
t=10;
Jvect=zeros(t,1);
for i=1:length(Jvect)
    [c_it, J_it]=predictor(rx, i);
    Jvect(i)=J_it;
end

% figure('Name', 'J over N');
% plot(1:t,Jvect);
% title('$J_{min}$ vs N');
% xlim([1 t]);
% xlabel('N'); ylabel('$J_{min}$');
%{
coeff=[1; copt];
A = tf([1 copt.'], 1,1);
figure, pzmap(A)
%}

% Variance of the AR model of the spectral lines
% upp_limit = 100;
% sigma_w = zeros(1, upp_limit);
% for N = 1:upp_limit
%    [~, sigma_w(N)] = findAR(N, rx);
% end
% figure, plot(1:upp_limit, 10*log10(sigma_w)), grid
% title('Variance of the AR model of the spectral lines')
% ylabel('$J_{min}$ [dB]'), xlabel('Order of the AR(N) model')
%% AR model
% Coefficients of Wiener filter
[a, s_white, d]=findAR(N, rx);
[H_w, omega] = freqz(1, [1; a], Nsamples, 'whole');

figure('Name','AR model estimate of the PSD');
plot(omega/(2*pi), 10*log10(s_white*abs(H_w).^2));
title('AR model estimate of the PSD');
xlabel('f');
ylabel('Amplitude (dB)');
ylim([-15 40]);

%% Final spectral plot
figure('Name', 'Spectral Analysis');
set(0,'defaultTextInterpreter','latex')          % to use LaTeX format
set(gca,'FontSize',10);
hold on;
plot((1:Nsamples)/Nsamples, 10*log10(abs(Pbt1)), 'Color', 'b')
plot((1:Nsamples)/Nsamples, 10*log10(Pper), 'g:')
plot((1:Nsamples)/Nsamples, 10*log10(Welch_P), 'r-.')
plot(omega/(2*pi), 10*log10(s_white*(abs(H_w)).^2), 'Color', 'm');
plot((1:Nsamples)/Nsamples, b, 'k:');
title('Spectral analysis');
legend('Correlogram', 'Periodogram', 'Welch', ['AR(' int2str(N) ')']);
hold off;
xlabel('Normalized frequency');
ylabel('Estimated PSD (dB)');
ylim([-15 30]);

[H, www] = freqz([1; a], 1, Nsamples, 'whole');
figure('Name', 'Z-plane for error predictor A(z)');
zplane([1; a].', 1);
title('Z-plane for error predictor A(z)');

%% Estimation of phases phi1 and phi2 based on coefficients of A(z)

% We notice that we have poles close to the unit circle at phases
% corresponding to phi2 and phi2
ph = angle(roots([1;a]))/(2*pi);