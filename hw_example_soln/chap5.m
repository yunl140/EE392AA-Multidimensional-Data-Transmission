clear;close all;clc;
%% 1+0.9D^-1
% H = [0.9 1 0; 0 0.9 1]/sqrt(.181);
% [~,~,M] = svd(H);
% Pm = M(:,1:2)*M(:,1:2)';
% Ctilde = Pm(:,1:2);
% Otilde=pinv(Ctilde)*Pm(:,3);
% Au_ut = [eye(2), Otilde];
% Rutut = Au_ut*Au_ut';
% [Q,~] = eig(Pm); % b/c Rxtxt = Pm
% Pq = Q(:,2:3)*Q(:,2:3)';
% A = Pq*Ctilde; % full rank, so u1 = utilde
% disp('Set Rvv = I:')
% Phi=sqrtm(Rutut);
% Ax_v = Phi\Au_ut;
% Ax_v*Ax_v' % check white v
% Htilde = H*pinv(Ax_v); % y = Htilde*v+n
% Rf = Htilde'*Htilde;
% Rbinv = Rf + eye(2);
% G = chol(Rbinv);
% S0 = diag(G).^2
% SNR = (prod(S0))^(1/3)-1 % matches
% 
% disp('Set Rvv to diagonal matrix but not I')
% Phi2 = [Phi(:,1), Phi(:,2)*2];
% Rvv = diag([1, .25]);
% Rutut - Phi2*Rvv*Phi2' % checks
% Rb2inv = Rf + inv(Rvv);
% Rb2 = inv(Rb2inv);
% Rzz = Rf*Phi2*Rvv*Phi2'*Rf+Rf;
% Ree = Rvv-Rb2*Rzz*Rb2; % does not match Rb2
% G2 = chol(Rb2inv);
% S02 = diag(G2).^2;

%% Example 5.1.4
% H2 = toeplitz([1 0], [1 0.9 0])/sqrt(.181);
% H1 = toeplitz([1 0], [1 -1 0])/sqrt(.181);
% disp('Rxx=I, order [2,1], process together:');
% H = [H2, H1];
% sum_rate = 0.5*log2(det(H*eye(6)*H' + eye(2))); % reverse order gives same rate
% fprintf('-sum rate: %.4f\n', sum_rate);
% Rbinv = H'*H + eye(6);
% [G,S0] = ldl(Rbinv, 'upper');
% ind_rate = 0.5*log2(diag(S0));
% fprintf('-individual rates: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n\n', ind_rate);
% forward_filter = inv(S0)*inv(G')*H'
% backward_filter = G
% 
% disp('optimize input Rxx,u individually, order [2,1], flat energy on used dimensions:');
% M2 = orth(H2');
% Rxx2 = 1.5*(M2*M2')
% H2t = H2*M2*sqrt(1.5);
% M1 = orth(H1');
% Rxx1 = 1.5*(M1*M1')
% H1t = H1*M1*sqrt(1.5);
% Ht = [H2t, H1t];
% Rxxt = blkdiag(Rxx2, Rxx1);
% sum_rate = 0.5*log2(det(Ht*eye(4)*Ht' + eye(2)));
% fprintf('-sum rate: %.4f\n', sum_rate);
% Rbinv = Ht'*Ht + eye(4);
% G = chol(Rbinv);
% S0 = diag(diag(G).^2);
% G = G./diag(G);
% ind_rate = 0.5*log2(diag(S0));
% fprintf('-individual rates: [%.4f, %.4f, %.4f, %.4f]\n\n', ind_rate);
% forward_filter = inv(S0)*inv(G')*Ht'
% backward_filter = G

%% Example 5.2.1
% H = toeplitz([.9 0], [.9 1 0])/sqrt(.181);
% M = orth(H');
% Pm = M*M';
% C = Pm(:,1:2);
% O = [eye(2), pinv(C)*Pm(:,3)];
% Ruu = O*O'
% G_phi = lohc(Ruu)
% Ht = H*C*G_phi
% Rbinv = Ht'*Ht + eye(2);
% [G,S0] = ldl(Rbinv, 'upper');
% snr = det(S0)^(1/3) - 1;
% fprintf('SNR_GDFE,U = %.4f dB\n', 10*log10(snr));
% forward_filter = inv(S0)*inv(G')*Ht'./sqrt(.181)
% backward_filter = G

%% Example 5.2.2
N = 8;
disp('---- CDFE ----')
H = toeplitz([.9, zeros(1,N-2), 1], [.9, 1, zeros(1,N-2)])/sqrt(.181);
Rbinv = H'*H+eye(N);
[G,S0] = ldl(Rbinv, 'upper');
snr = exp(sum(log(diag(S0)))/(N+1)) - 1;
fprintf('SNR_CDFE,U = %.4f dB \n', 10*log10(snr));
forward_filter = inv(S0)*inv(G')*H'./sqrt(.181);
backward_filter = G;
disp('---- GDFE w/ increased flat energy ----')
H = toeplitz([.9, zeros(1, N-1)], [.9 1 zeros(1, N-1)])/sqrt(.181)*sqrt((N+1)/N);
M = orth(H');
Pm = M*M';
C = Pm(:,1:N);
O = [eye(N),pinv(C)*Pm(:,end)];
Ruu = O*O';
G_phi = lohc(Ruu);
Ht = H*C*G_phi;
Rbinv = Ht'*Ht+eye(N);
[G,S0] = ldl(Rbinv, 'upper');
snr = exp(sum(log(diag(S0)))/(N+1)) - 1;
fprintf('SNR_GDFE,U = %.4f dB \n', 10*log10(snr));
