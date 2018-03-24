function [a, sw, det_R]=findAR(N, rx)
%finds the parameters of an AR model of a r.p.
%INPUT: order of the filter N, autocorrelation of the r.p. rx
%OUTPUT: coefficients of the filter in the vector a, value of the variance
%and the determinant of R det_R to check for ill-conditioning
%of the white process

col=rx(1:N);
row=conj(col);
%construct the autocorrelation matrix NxN
R=toeplitz(col, row);
%compute the vector r
r=rx(2:N+1);
%compute the coefficients
det_R=det(R);
a=-inv(R)*r;
sw=abs(R(1,1)+r'*a);
end