%% problem 2.30 (e)
%clear; clc;
% noise-whitened channel
H = [8 2 6;4 2 3;2 2 1.5]*10;
Rxx = diag([1,4,1]);
%H = [8 6 4;6 4.5 3;2 2 2]*10;
%Rxx = diag([3 4 2]);
%[Rwcn, sumrate] = wcnoise(Rxx, H, 1);
%Swcn = inv(Rwcn) - inv(H*H'+ Rwcn) % -> confirm 2nd user is secondary

% reorder users so that user 3 is secondary
H = H([1,3,2],[1,3,2]);
Rxx = diag([1,1,4]);
U0 = 2;

reverse_idx_func = @(X) X(end:-1:1, end:-1:1);
% temporarily remove secondary user
H1 = H(1:2,:);
[Rwcn, b] = wcnoise(Rxx, H1, 1, 1e-5);
Swcn = inv(Rwcn) - inv(H1*Rxx*H1'+ Rwcn);
[~,Swcn1,Qwcn] = svd(Swcn);
tmp = Qwcn/Rwcn*H1;
%tmp=H1;
[Q,R] = qr(reverse_idx_func(tmp)');
Q = Q(end:-1:1, U0:-1:1);
R = R(U0:-1:1, U0:-1:1)';
Phi = reverse_idx_func(chol(reverse_idx_func(Q'*Rxx*Q)))';
A = Q*Phi;
Rb = inv(eye(2) + A'*H1'/Rwcn*H1*A);
DA = diag(diag(R*Phi));
G = DA\R*Phi;
S0 = DA/Swcn1*DA;
rates = log2(diag(S0));
sum(rates)

% energize secondary user
E_p = 1/5;
E_s = 4/5;
rates_p = log2(diag(S0*E_p));
% channel for the secondary receiver w.r.t. v3
Htilde = H*A;
H3 = Htilde(3,:);
snr_s = sum(H3)^2*E_s/(1+norm(H3)^2*E_p);
rate_s = log2(1 + snr_s)
sum(rates_p) + rate_s

