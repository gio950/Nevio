close all; clear global; clearvars; clc;

%% GENERATE THE PROCESS x(k), 1 REALIZATION

Nsamples=800;
% Frequencies of the exponentials
f1=0.17;
f2=0.78;
% Generate the white noise (2 components)
sigmaw=0.1;
x = zeros(800,300);

for i=1:300
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
    x(:,i) = xi + 1i*xq;
end
save('realizations.mat', 'x');