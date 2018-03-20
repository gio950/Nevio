function [rx]=autocorrelation_Unb(x)
% Unbiased autocorrelation estimator
%INPUT: r.p. x, length of the autocorrelation Lcorr
%OUTPUT: autocorrelation estimate vector rx of length K=length(x)
%every index is augmented by 1 because matlab starts from 1 and not 0
K=length(x);

rx=zeros(K, 1);
   for n=1:K
      %first x that has k as argument
      xnk=x(n:K);
      %second x that has k-n as argument and the conjugate
      xconj=conj(x(1:(K-n+1)));
      rx(n)=(xnk.'*xconj)/(K-n+1);
   end
end