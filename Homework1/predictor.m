function [copt, Jmin, det_R]=predictor(rx, N)

[a, sw ]=findAR(N, rx);
copt=-a;
Jmin=sw;

end