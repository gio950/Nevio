function [ autoc_full ] = full_autocorr(x, rx)
% Finds the negative symmetrical part of the autocorrelation estimate 
%INPUT: r.p. x, autocorrelation estimate rx
%OUTPUT: vector autoc_full with the autocorrelation (both negative and positive terms) 

L=length(rx);
K = length(x);
% make autocorrelation simmetric
autoc_full = zeros(K, 1);
autoc_full(1:L) = rx;
temp = flipud(conj(rx));
% it's the same as puttig the conjugate, flipped, at the end of this
% vector, since fft is periodic
autoc_full((K-L+1):K) = temp(1:length(temp));
end