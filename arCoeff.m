function [ c ] = arCoeff(autocorr, N)
% autocorr: estimated autocorrelation of input signal
col=autocorr(1:N);
row=conj(col);
R=toeplitz(col,row);
% vector r: one sample forward with respect to matrix R
r = autocorr(2:N+1);

c = inv(R) * r;
end