clear;clc;close all;
%% P1
disp('-------- 1.b --------')
fprintf('one-sided PSD level: -60 dBm/Hz\n\n');
disp('-------- 1.c --------')
ra = 8*20*log2(1+200);
fprintf('Design A best rate: %.4f Mbps\n\n', ra);
disp('-------- 1.d --------')
rb = 2*20*log2(1+2*200);
fprintf('Design B best rate: %.4f Mbps\n\n', rb);
disp('-------- 1.e --------')
ra = 8*20*log2(1+200/10^.3);
fprintf('Design A rate w/ 3dB gap: %.4f\n', ra);
snr_geo = sqrt(2*200+1)-1;
rb = 2*40*log2(1+snr_geo/10^.3);
fprintf('Design B rate w/ 3dB gap: %.4f\n\n', rb);
disp('-------- 1.f --------')
snr_geo = sqrt(200/2+1)-1;
rc = 320*log2(1+snr_geo/10^.3);
fprintf('Design Proud rate w/ 3dB gap: %.4f\n\n', rc);
%% P2
h = [0.9 1];
N0 = 2*.181;
Es = [4 32];
htilde = h/sqrt(N0).*sqrt(Es);
bbars = [1.5 1];
disp('-------- 2.a --------')
Gammas_top = htilde.^2./(2.^(2*bbars) - 1);
Gammas_bottom = htilde.^2./(1+htilde([2,1]).^2)./(2.^(2*bbars) - 1)
% produces [0.0219    2.0665]
% user 2 has nagative gap if it's at the bottom of the order
% so the order must be [2,1]
Gammas_db = 10*log10([Gammas_top(1), Gammas_bottom(2)]);
[Gammas_top(1), Gammas_bottom(2)]
fprintf('Gap is %.4f dB for user 2; %.4f dB for user 1\n\n', Gammas_db)
disp('-------- 2.b --------')
fprintf('<P> = 1e-6\n\n')
disp('-------- 2.d --------')
b2max = log2(1+htilde(1)^2);
b1max = log2(1+htilde(2)^2/(1+htilde(1)^2));
fprintf('max achievable rates: b2=%.4f, b1=%.4f\n\n', b2max, b1max);
disp('-------- 2.e --------')
bmax = log2(1+36/N0);
fprintf('max sum rate: %.4f bits/subsymbol\n\n', bmax);
disp('-------- 2.h --------')
% basically, SNR doubled
N0 = .181;
htilde = h/sqrt(N0).*sqrt(Es);
% order [2,1]
b2max = log2(1+htilde(1)^2);
b1max = log2(1+htilde(2)^2/(1+htilde(1)^2));
fprintf('order [2,1]: b2=%.4f, b1=%.4f, sum=%.4f\n', b2max, b1max, b2max+b1max);
% order [1,2]
b1max = log2(1+htilde(2)^2);
b2max = log2(1+htilde(1)^2/(1+htilde(2)^2));
fprintf('order [1,2]: b2=%.4f, b1=%.4f, sum=%.4f\n\n', b2max, b1max, b2max+b1max);
disp('-------- 2.i --------')
N0 = .181/4;
htilde = h/sqrt(N0).*sqrt(Es);
% order [2,1]
b2max = log2(1+htilde(1)^2);
b1max = log2(1+htilde(2)^2/(1+htilde(1)^2));
fprintf('order [2,1]: b2=%.4f, b1=%.4f, sum=%.4f\n', b2max, b1max, b2max+b1max);
% order [1,2]
b1max = log2(1+htilde(2)^2);
b2max = log2(1+htilde(1)^2/(1+htilde(2)^2));
fprintf('order [1,2]: b2=%.4f, b1=%.4f, sum=%.4f\n\n', b2max, b1max, b2max+b1max);
