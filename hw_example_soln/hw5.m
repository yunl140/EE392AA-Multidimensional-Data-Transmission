clear;close all;clc;
%% 2.29
disp('-------- 2.29 (c) --------')
disp('see figure 1')
cb1 = @(a) .5*log2(1+800*a);
cb2 = @(a) .5*log2(1+(1-a)./(.005+a));
figure(1)
fplot(@(a) cb1(a), @(a) cb2(a), [0,1]);

%% 2.30
clear;
disp('-------- 2.30 (a) --------')
Htilde = 10*[8 2 6;4 2 3;2 2 1.5]
disp('-------- 2.30 (b) --------')
Rxx = eye(3);
[Rwcn, bbar_sum] = wcnoise(Rxx, Htilde, 1)
disp('Verify primary/secondary users:')
Swcn = inv(Rwcn) - inv(Htilde*Rxx*Htilde' + Rwcn)
disp('-------- 2.30 (c) --------')
disp('Change Rxx to Rxx = ')
tmp = rand(3);
Rxx = tmp*tmp';
Rxx = 3*Rxx/trace(Rxx)
[Rwcn, bbar_sum] = wcnoise(Rxx, Htilde, 1)
disp('Verify primary/secondary users:')
Swcn = inv(Rwcn) - inv(Htilde*Rxx*Htilde' + Rwcn)
disp('-------- 2.30 (d) --------')
% reorder users
Hd = Htilde([1 3 2], [1 3 2]);
Rxx = diag([1 1 4]);
% temporarily remove secondary user
Hdprime = Hd(1:2,:);
[Rwcn, bbar_sum] = wcnoise(Rxx, Hdprime, 1, 1e-5);
[~,Swcn1, Qwcn] = svd(inv(Rwcn) - inv(Hdprime*Rxx*Hdprime' + Rwcn));
[R,Q] = rq(Qwcn/Rwcn*Hdprime, 0);
Phi = lohc(Q'*Rxx*Q);
A = Q*Phi;
DA = diag(R*Phi);
S0 = DA.^2./diag(Swcn1);
rates = log2(S0);
fprintf('If user 3 is not energized, r = [%.4f %.4f 0], sum rate: %.4f\n', rates, sum(rates));
% energize user 3
E_p = 1/5;
E_s = 4/5;
rates_p = log2(1+(S0-1)*E_p);
H3 = Hd*A;
H3 = H3(3,:);
rate_s = log2(1+sum(H3)^2*E_s/(1+H3*H3'*E_p));
fprintf('If user 3 is energized, r = [%.4f %.4f %.4f], sum rate: %.4f\n\n', rates_p, rate_s, sum(rates_p)+rate_s);

%% 2.31
clear;
disp('-------- 2.31 (b) --------')
H = [8 2 6;4 2 3; 2 1 1.5]*10;
H = H([2 3 1], [2 3 1])
H22 = H(1:2, 1:2);
[u,s,v] = svd(H22); % xtilde = v(1,:)*x, now 1D, y = u(:,1)*s*xtilde
Htilde = [H(:,1:2)*v(1,:)', H(:,3)]
disp('-------- 2.31 (d) --------')
% order [2 1]
rate_top2 = log2(1+svd(Htilde([1,2],1))^2);
rate_bottom1 = log2(1+Htilde(3,2)^2/(1+Htilde(3,1)^2));
fprintf('order [2 1]: r2= %.4f, r1= %.4f, sum= %.4f\n', rate_top2, rate_bottom1, rate_top2+rate_bottom1);
% order [1 2]
rate_top1 = log2(1+Htilde(3,2)^2);
rate_bottom2 = log2(1+svd(Htilde([1,2],1))^2/det(eye(2)+Htilde(1:2,2)*Htilde(1:2,2)'));
fprintf('order [1 2]: r2= %.4f, r1= %.4f, sum= %.4f\n\n', rate_bottom2, rate_top1, rate_bottom2+rate_top1);

%% 2.32
clear;
g1upper = 4/2; g1lower = 3/2;
g2upper = 3/2; g2lower = 4/2;
% order [2,1] is always superior to order [1,2] in the BC part
rates_BC = @(x) 0.5*[log2(1+g1upper*x); log2(1+g1lower*(1-x)./(1+g1lower*x))];
x = 0:0.001:1;
rbc = rates_BC(x);
% now the MAC pentagone
rmac = 0.5*[log2(1+g2upper), log2(1+g2lower), ...
    log2(1+g2upper/(1+g2lower)), log2(1+g2lower/(1+g2upper))];
figure(2); hold on;
plot(rbc(1,:), rbc(2,:));
plot([rmac(1), rmac(1), rmac(3), 0], [0, rmac(4), rmac(2), rmac(2)])
legend('BC rates/dim', 'MAC rates/dim');
xlabel('b-upper');
ylabel('b-lower');
box on;
grid on;
hold off;