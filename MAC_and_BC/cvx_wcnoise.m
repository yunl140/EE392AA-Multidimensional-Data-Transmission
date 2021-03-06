function [Rnn, sumRatebar] = cvx_wcnoise(Rxx, H, Lyu)
%cvx_wcnoise This function computes the worst-case noise for a given input
%autocorrelation Rxx and channel matrix.
%   Arguments:
%       - Rxx: input autocorrelation, size(Lx, Lx)
%       - H: channel response, size (Ly, Lx)
%       - Lyu: number of antennas at each user, scalar/vector of length U
%   Outputs:
%       - Rnn: worst-case noise autocorrelation, with white local noise
%       - sumRatebar: maximum sum rate/real-dimension

eps_margin = 1e-9;
[Ly,Lx] = size(H);
if length(Lyu) == 1
    U = Ly/Lyu;
    Lyus = ones(1,U)*Lyu;
else
    U = length(Lyu);
    Lyus = Lyu;
end
cum_Lyu = cumsum(Lyus);
cum_Lyu = [0, cum_Lyu(1:end-1)];
us = 1:U;

[Q,D,~]=svd(Rxx);
rD = rank(D);
sqD = sqrt(diag(D));
sqD(rD+1:end) = 0;
Htilde = H*Q*diag(sqD);
cvx_begin sdp quiet
cvx_solver mosek
cvx_precision high
variable Rnn(Ly, Ly) hermitian semidefinite
variable Z(Lx, Lx) hermitian semidefinite
dual variable S1{U}
dual variable S2{U}
maximize log_det(eye(Lx) - Z)
subject to
for u = us
    S1{u}: eye(Lyus(u))==real(Rnn(cum_Lyu(u) + (1:Lyus(u)), cum_Lyu(u) + (1:Lyus(u))));
    S2{u}: 0==imag(Rnn(cum_Lyu(u) + (1:Lyus(u)), cum_Lyu(u) + (1:Lyus(u))));
end
Rnn - eps_margin*eye(Ly) ==  hermitian_semidefinite(Ly);
[Rnn + Htilde*Htilde', Htilde; Htilde', Z] == hermitian_semidefinite(Ly+Lx);
cvx_end

sumRatebar = -0.5*log2(exp(1))*cvx_optval;
