%% 2.26
clear;close all;clc;
disp('------------ 2.26 ------------')
% part b
disp('------ part b ------')
u = [sqrt(.75), .5; -.5, sqrt(.75)];
Rnn = u*diag([9/256, 1/16])*u'
sqrtRnn = u*diag([3/16, 1/4])
sqrtRnn = sqrtm(Rnn)
H = [5 2 1;3 1 1];
Htilde = inv(sqrtRnn)*H
% part c
disp('------ part c ------')
Ly=2;U=3;
cvx_begin
    variable x(U)
    maximize .5*log2(exp(1))*log_det(eye(Ly) + Htilde*diag(x)*Htilde')
    subject to
    x <= 2;
    x >= 0;
cvx_end
% part d
disp('------ part d ------')
0.5*log2(det(Htilde*2*eye(3)*Htilde' + eye(2)))
% part e
disp('------ part e ------')
rates2 = @(h) 0.5*log2(diag(chol(h'*h+eye(2))).^2);
Htilde = Htilde*sqrt(2);
fprintf('E3=0, order [2,1], rates: %.4f %.4f, sum: %.4f\n', rates2(Htilde(:,2:3)), sum(rates2(Htilde(:,2:3))));
fprintf('E3=0, order [1,2], rates: %.4f %.4f, sum: %.4f\n', rates2(Htilde(:,3:-1:2)), sum(rates2(Htilde(:,3:-1:2))));
fprintf('E2=0, order [3,1], rates: %.4f %.4f, sum: %.4f\n', rates2(Htilde(:,[1,3])), sum(rates2(Htilde(:,[1,3]))));
fprintf('E2=0, order [1,3], rates: %.4f %.4f, sum: %.4f\n', rates2(Htilde(:,[3,1])), sum(rates2(Htilde(:,[3,1]))));
fprintf('E1=0, order [3,2], rates: %.4f %.4f, sum: %.4f\n', rates2(Htilde(:,1:2)), sum(rates2(Htilde(:,1:2))));
fprintf('E1=0, order [3,2], rates: %.4f %.4f, sum: %.4f\n\n', rates2(Htilde(:,[2,1])), sum(rates2(Htilde(:,[2,1]))));
disp('------ part f ------')
rates3 = @(H, order) 0.5*log2(diag(chol(H(:,order)'*H(:,order)+eye(3))).^2);
vertices = zeros(6,3);
idx = 1;
for order = perms(1:3)'
    r = rates3(Htilde,order');
    r = r(order');
    vertices(idx,:) = r;
    fprintf('order [%1i %1i %1i], rates %.4f %.4f %.4f, sum: %.4f\n', order', r, sum(r));
    idx = idx + 1;
end
disp('------ part g ------')
inhull([.1, 2.4, 3.2], vertices)
disp('------ part h ------')
Htilde = inv(sqrtRnn)*H;
Ly=2;U=3;
cvx_begin
    variable x(U)
    maximize .5*log2(exp(1))*log_det(eye(Ly) + Htilde*diag(x)*Htilde')
    subject to
    sum(x) <= 6;
    x >= 0;
cvx_end

%% 2.27
invH1general = @(f,a,s) s^2./(1+a^2 - 2*a*cos(2*pi*f));
invH2general = @(f,a,s) s^2./(1+a^2 + 2*a*cos(2*pi*f));
invH1 = @(f) .181./(1.81 - 1.8*cos(2*pi*f));
invH2 = @(f) .181./(1.81 + 1.8*cos(2*pi*f));
figure(1); hold on;
fplot(invH1, [0,0.5]);
fplot(invH2, [0,0.5]);
legend('sigma^2/|H_1(f)|^2', 'sigma^2/|H_2(f)|^2');
xlabel('f');
grid on;
box on;
hold off;

figure(2); hold on;
h1=fplot(@(f) invH1general(f, 0.5, .181), [0,0.5]);
h2=fplot(@(f) invH2general(f, 0.5, .181), [0,0.5]);
legend([h1,h2],'sigma^2/|H_1(f)|^2', 'sigma^2/|H_2(f)|^2');
fs = 0.25:0.01:0.5;
lambda = 0.1;
UB = lambda*ones(1,length(fs));
x2 = [fs, fliplr(fs)];
inBetween = [UB, fliplr(invH1general(fs, 0.5, .181))];
fill(x2, inBetween, 'b');
fs = 0:0.01:0.25;
x2 = [fs, fliplr(fs)];
inBetween = [UB, fliplr(invH2general(fs, 0.5, .181))];
fill(x2, inBetween, 'r');
xlabel('f');
grid on;
box on;
hold off;
lambda = 2 + 4*.181/pi/sqrt(1.81^2-1.8^2)*atan(1/19)
rate = integral(@(f) log2(lambda/.181*(1.81+1.8*cos(2*pi*f))), 0, 0.25)
sum_rate = rate*2