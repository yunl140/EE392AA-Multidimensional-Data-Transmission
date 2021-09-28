clear;clc;close all;
%% P5.5
Rxx = [1 .5 1;.5 1 .5; 1 .5 1];
% part a
disp('-------- 1(a) --------')
fprintf('N* = %d, N* = Nbar*\n', rank(Rxx));
disp('-------- 1(b) --------')
[V,D] = eig(Rxx)
Q = V(:,2:3)
Pq = Q*Q'
disp('-------- 1(c) --------')
disp('solution with svd:')
[u,s,v]=svd(Rxx);
A = u(:,1:2)*sqrtm(s(1:2,1:2))
Ruu = pinv(A)*Rxx*pinv(A)'
disp('solution with Qtilde:')
A=Q
Ruu = pinv(A)*Rxx*pinv(A)'

%% P5.8
clear;
disp('-------- 3(a) --------')
H = 1/sqrt(.181) * toeplitz([.9 0], [.9 1 0]);
[~,~,M] = svd(H);
Mtilde = -M(:,1:2);
Pm = Mtilde*Mtilde'; % xt = Pm*x
Ctilde = Pm(:,1:2);
Otilde = pinv(Ctilde)*Pm(:,3);
Oextend = [eye(2), Otilde]; % ut = Oextend*u = Oextend*x
Rutut = Oextend*Oextend'
disp('-------- 3(b) --------')
us = 1 - 2*(dec2bin(0:7)'-'0')
utildes = Oextend*us
disp('-------- 3(c) --------')
Rxtxt = Pm;
fprintf('wasted energy: %.f / 3 units\n', 3-trace(Rxtxt));
disp('-------- 3(d) --------')
xs = Ctilde*utildes
disp('-------- 3(e) --------')
y = H*us

%% P5.9
clear;
disp('-------- 4(a) --------')
N = 3; v = 1;
H = toeplitz([1 zeros(1,N-2) 1], [1 1 zeros(1,N-2)])/sqrt(.001);
Rf = H'*H;
Rbinv = Rf+eye(N);
[sqrtRxx,S0] = ldl(Rbinv, 'upper');
SNR_CDFEu = det(S0)^(1/(N+v))-1;
W = inv(S0)*inv(sqrtRxx');
Pe = 2*(15/16)*qfunc(sqrt(3*SNR_CDFEu/(16^2-1)));
sqrtRxx
W
fprintf('SNR_CDFE,U = %.4f = %.4f dB, Pe = %.4f\n\n',...
    SNR_CDFEu, 10*log10(SNR_CDFEu), Pe)

disp('-------- 4(b) --------')
[SNR_DFEu,~,~]=dfecolor(1, [1 1], 3, 1, -1, 1, 0.001);
Pe = 2*(7/8)*qfunc(sqrt(3*10^(SNR_DFEu/10)/(8^2-1)))
Peprop = Pe*8

disp('-------- 4(c) --------')
H = toeplitz([1 0 0], [1 1 0 0])/sqrt(.001);
[~,~,M] = svd(H);
Ctilde = M(:,1:3);
Htilde = H*Ctilde*sqrt(4/3);
Rbinv = Htilde'*Htilde+eye(3);
[sqrtRxx,S0] = ldl(Rbinv, 'upper');
SNR_CDFEu = det(S0)^(1/4)-1;
W = inv(S0)*inv(sqrtRxx');
Pe = 2*7/8*qfunc(sqrt(3*SNR_CDFEu/(16^2-1)));
sqrtRxx
W
fprintf('SNR_CDFE,U = %.4f = %.4f dB, Pe = %.4f\n\n',...
    SNR_CDFEu, 10*log10(SNR_CDFEu), Pe);

disp('-------- 4(d) --------')
Ns = [10, 20, 100, 200];
for N = Ns
    H = toeplitz([1 zeros(1, N-1)], [1 1 zeros(1, N-1)])/sqrt(.001);
    [~,~,M] = svd(H);
    Htilde = H*M(:,1:end-1)*sqrt((N+1)/N);
    S0 = diag(Htilde'*Htilde) + 1;
    snr = exp(sum(log(S0))/(N+1)) - 1;
    m = floor(sqrt(3*snr/(qfuncinv(1e-6/2)^2) + 1));
    fprintf('N = %d, SNR = %.4f = %.4f, %d PAM\n', N, snr, 10*log10(snr), m);
end

%% P5.11
disp('-------- 5(a) --------')
H = [1 .5 .5; .5 1 .5; .5 .5 1]*10;
[V,gs] = eig(H'*H);
gs = diag(gs);
Ex = (3 + sum(1./gs))/3 - 1./gs;

Rxx = V*diag(Ex)*V'
[Gx,Sx]=ldl_rev_idx(Rxx)
disp('-------- 5(b) --------')
sqrtRxx = lohc(Rxx);
Htilde = H*sqrtRxx;
Rf = Htilde'*Htilde
Rbinv = Rf + eye(3)
disp('-------- 5(c) --------')
[G, S0] = ldl(Rbinv, 'upper');
G
W = inv(S0)*inv(G')
disp('-------- 5(d) --------')
snr = det(S0)^(1/3)-1;
snr_db = 10*log10(snr);
bbar = 0.5*log2(1+snr);
fprintf('SNR = %.4f dB, bbar = %.4f bits/dimension\n', snr_db, bbar);
