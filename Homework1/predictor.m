function [copt, Jmin, det_R]=predictor(rx, N)

[a, sw, det_R ]=findAR(N, rx);
copt=-a;
Jmin=sw;

end