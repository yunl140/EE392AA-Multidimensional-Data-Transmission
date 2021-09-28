clear; clc; close all;
%% 5.16
% % (a)
% disp('-------- 5.16 (a) --------')
h = cat(3, [1 -.5; .9 1], [1 -.4; 0 -1])*10;
H = fft(h, 8, 3)/sqrt(8)
% % (b)
% disp('-------- 5.16 (b) --------')
% Hcell = num2cell(H, [1,2]);
% gains = reshape(cell2mat(cellfun(@(h) svd(h).^2, Hcell, 'UniformOutput', 0)), 1, 16);
% % waterfill
% gains = sort(gains, 'descend');
% for Gstar = 16:-1:1
%     K = (2*(8/9) + sum(1./gains(1:Gstar)))/Gstar;
%     if K - 1/gains(Gstar)>0
%         break
%     end
% end
% Es = zeros(1,16);
% Es(1:Gstar) = K - 1./gains(1:Gstar)
% sum_rate = sum(log2(K*gains(1:Gstar)));
% fprintf('maximum data rate is %.4f bits/DMT-symbol\n', sum_rate);
% % (c)
% disp('-------- 5.16 (c) --------')
% bs = zeros(2,8);
% H = H;
% for n = 1:8
%     Hn = H(:,:,n)*sqrt(1/9);
%     Rbinv = Hn'*Hn + eye(2);
%     [G,S0,P] = ldl(Rbinv, 'upper');
%     b(:,n) = log2(diag(S0));
% end
% b
% sum_rate_mac = sum(b,'all');
% fprintf(' user 2: %.4f bits/symbol\n user 1: %.4f bits/symbol\n sum: %.4f bits/symbol\n',...
%     sum(b,2), sum_rate_mac);
% % (d)
% disp('-------- 5.16 (d) --------')
% Hcell = num2cell(H,[1]);
% Hexpand = [blkdiag(Hcell{1,1,:}), blkdiag(Hcell{1,2,:})];
% cvx_begin
%     variable Es(16)
%     maximize log2(exp(1))*log_det(eye(16) + Hexpand*diag(Es)*Hexpand')
%     subject to
%         sum(Es) <= 16/9;
%         Es >= 0;
% cvx_end
% sum_rate_esum_mac = cvx_optval
% gamma_mac = (2^(sum_rate/9)-1)/(2^(sum_rate_esum_mac/9)-1);
% fprintf('gamma_mac = %.4f = %.4f dB\n\n', gamma_mac, 10*log10(gamma_mac));

%% 5.17
% % (a)
% disp('-------- 5.17 (a) --------')
% D=exp(-1i*(2*pi/8)*(0:7));
% H22=ones(1,8)+D;
% H21=-.5*ones(1,8)-.4*D;
% H12=.9*ones(1,8);
% H11=ones(1,8)-D;
% 
% scale=sqrt(100/8);
% 
% H0=scale*[H22(1) H21(1)
% H12(1) H11(1)] ;
% H1=scale*[H22(2) H21(2)
% H12(2) H11(2)];
% H2=scale*[H22(3) H21(3)
% H12(3) H11(3)] ; 
%  H3=scale*[H22(4) H21(4)
% H12(4) H11(4)] ; 
%  H4=scale*[H22(5) H21(5)
% H12(5) H11(5)] ;
% H5=scale*[H22(6) H21(6)
% H12(6) H11(6)] ;
% H6=scale*[H22(7) H21(7)
% H12(7) H11(7)] ;
% H7=scale*[H22(8) H21(8)
% H12(8) H11(8)]  ;
% 
% H=zeros(2,2,8);
% H(:,:,1)=H0;
% H(:,:,2)=H1;
% H(:,:,3)=H2;
% H(:,:,4)=H3;
% H(:,:,5)=H4;
% H(:,:,6)=H5;
% H(:,:,7)=H6;
% H(:,:,8)=H7;
% 
% [E, b, theta, w, ~] = admMAC(H, [2; 1], 8/9*[1; 1]);
% E
% sum(E,2)
% b=real(b)
% sum(b,2)
% % (b)
% disp('-------- 5.17 (b) --------')
% [E, b, theta, w, ~] = admMAC(H, [16; 8], 8/9*[1; 1])
% [E, theta, b] = minPMAC(H, [16; 8], [1; 1])

%% 5.18
% % (a)
% disp('-------- 5.18 (a) --------')
% Hinit = [toeplitz([1 0], [1 1 0]), toeplitz([-.5 0], [-.5 -.4 0]);
%     toeplitz([.9 0], [.9 0 0]), toeplitz([1 0], [1 -1 0])]
% % (b)
% disp('-------- 5.18 (b) --------')
% [Ly,Lx] = size(Hinit);
% L = lcm(Ly,Lx)
% % (c)
% disp('-------- 5.18 (c) --------')
% %Htilde = [Hinit,Hinit;Hinit,Hinit;Hinit,Hinit]
% Hinit = Hinit*10;
% Rxxmac = repmat(eye(3)/3, 1, 1, 2);
% Hmac = cat(3, Hinit(:,1:(Lx/2)), Hinit(:,(Lx/2)+1:end));
% Rxxbc_u = mac2BcMimo(Rxxmac, Hmac, [2;1])
% % (d)
% disp('-------- 5.18 (d) --------')
% Rxxbc = sum(Rxxbc_u,3)

%% 5.19
% % (a)
% disp('-------- 5.19 (a) --------')
% reverse_idx = @(X) X(end:-1:1, end:-1:1,:);
% Hbc = conj(reverse_idx(permute(H,[2 1 3])))
% 
% % (b)
% disp('-------- 5.19 (b) --------')
% % align all indices conventions
% Hbc = conj(permute(H,[2,1,3]));
% 
% Rxxbaru = zeros(2,2,2,8);
% Rxxbar = zeros(2,2,8);
% bmac = zeros(2,8);
% bbc = zeros(2,8);
% order = [1,2];
% for n=1:8
%     Hn = reshape(H(:,:,n),2,1,2);
%     S = ones(1,1,2)/9;
%     Rxxbaru(:,:,:,n) = mac2BcMimo(S, Hn, order);
%     Rxxbar(:,:,n) = sum(Rxxbaru(:,:,:,n),3);
%     % mac rates
%     bmac(:,n) = MACrates(H(:,:,n),eye(2)/9,order);
%     % bc rates
%     bbc(:,n) = BCrates(Hbc(:,:,n), Rxxbaru(:,:,:,n), order);
% end
% 
% bmac
% bbc = real(bbc)

%% 5.20
% % (a)
% disp('-------- 5.20 (a) --------')
% G = [50 50; 0 1].^2;
% Es = [1 1];
% eps = 1;
% while eps > 1e-5
%     g1 = G(1,1)/(1+G(1,2)*Es(2));
%     g2 = G(2,2)/(1+G(2,1)*Es(1));
%     [gs, idx]=sort([g1,g2],'descend');
%     for Nuse = 2:-1:1
%         K = (2 + sum(1./gs(1:Nuse)))/Nuse;
%         Es_new = K-1./gs(1:Nuse);
%         if Es_new(end) > 0
%             break
%         end
%     end
%     Es_new(idx) = Es_new;
%     eps = max(abs(Es_new - Es));
%     Es = Es_new;
% end
% bs = log2(1+Es.*[g1,g2]);
% fprintf('E1 = %.4f, b1 = %.4f\nE2 = %.4f, b2 = %.4f\nrate sum = %.4f\n\n',...
%     Es(1), bs(1), Es(2), bs(2), sum(bs));
% 
% % (b)
% disp('-------- 5.20 (b) --------')
% % note that user 2 is intereference-free
% % order [2,1] at user 1
% g = G(1,1);
% E1 = 1 - 1/g; E2 = 2 - E1;
% bs_21 = [log2(1+g*E1), log2(1+E2)];
% bsum_21 = sum(bs_21)
% % order [1,2] at user 1
% %E2 = (sqrt(2+2*g^2) - 1)/g;
% E2 = 0;
% E1 = 2 - E2;
% 
% bs_12 = [log2(1 + g*E1/(1+g*E2)), log2(1+E2)];
% bsum_12 = sum(bs_12)
% fplot(@(E2) log2(1+(2-E2)./(1/g+E2)) + log2(1+E2), [0,2])
% figure(2)
% fplot(@(E1) log2(1+2500*E1)+log2(2-E1), [0, 2499/2500]); hold on;
% fplot(@(E1) log2(2500*(2-E1)), [2499/2500,2]); hold off;
% figure(3)
% fplot(@(E1) log2(1+2-E1), [0,2]); hold on;
% fplot(@(E1) log2(1+(2-E1)./(1/g+E1)), [0,2]); hold off;
% legend('b2','b2 ub')