function [corr] = correlogram( x, window, rx, L )
% Compute the PSD estimate using the correlogram method
%INPUT: r.p. x, window (hamming or rect, passed as a vector of samples
%starting from the origin), autocorrelation estimate rx, L number of
%samples used for the autocorrelation
%OUTPUT: corr, the PSD estimate using the correlogram method

K = length(x);
autoc_complete = full_autocorr(x, rx);
%computes the symmetric part of the window (negative samples)
full_win = zeros(K, 1);
%exploits periodicity instead of going to "negative indexes"
full_win(1 : L + 1) = window(L + 1 : 2*L + 1);
full_win(K - L + 1 : K) = window(1 : L);

windowed_autoc = autoc_complete .* full_win;
corr = fft(windowed_autoc);

end