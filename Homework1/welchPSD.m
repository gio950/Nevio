function [welch_est, Ns] = welchPSD(inputsig, window, overlaps)
% Length of the window
D = length(window);
% Length of input signal
K = length(inputsig);
% Normalized energy of the window
Mw = sum(window .^ 2) * (1/D);
% Number of subsequences
N_s = floor((K-D)/(D-overlaps) + 1);
%Initialization of each periodogram
P_per = zeros(K, N_s);

for s = 0:(N_s-1)
    % Windowed input
    x_s = window .* inputsig(s*(D-overlaps)+1:s*(D-overlaps)+D);
    % DFT on K samples of windowed input
    X_s = fft(x_s, K);
    % Periodogram for the window
    P_per(:,s+1) = (abs(X_s)).^2 * (1/(D*Mw)); % Tc = 1;
end
% Sum of all periodograms
welch_est = sum(P_per, 2) * (1/N_s);
Ns = length(welch_est);
end