R=@(sinr) 0.5*log2(1+sinr);
g11 = 240; g12 = 225; g21 = 100; g22 = 90;
weight2=1.05;
%% individual power constraint (1,1), ordering
% max sum rate achieved when p1=p2=1
% dec 1 then 2
R12 = R( min(g11/(1+g12), g21/(1+g22)) ) + weight2*R(g22);
% dec 2 then 1
R21 = R(g11) + weight2*R( min(g22/(1+g21), g12/(1+g11)) );

if R12 > R21
    fprintf('Individual power constraint, NOMA no time share:\nDecode 1 then 2\nUse power vector (%.4f, %.4f)\nSum Rate: %.4f\n\n', [p12, R12]);
else
    fprintf('Individual power constraint, NOMA no time share:\nDecode 2 then 1\nUse power vector (%.4f, %.4f)\nSum Rate: %.4f\n\n', [p21, R21]);
end
%% sum power constraint (1+1), ordering
% dec 2 then 1
[p1, fval] = fminbnd(@(p) -R(g11*p)-weight2*R(g12*(2-p)./(g11*p+1)), 0, 2);
fprintf('Sum power constraint, decode 2 then 1:\nPower: (%.4f, %.4f)\nSum Rate: %.4f\n', [p1, 2-p1, -fval]);
% dec 1 then 2
[p2, fval] = fminbnd(@(p) -weight2*R(g22*p)-R(g21*(2-p)./(g22*p+1)), 0, 2);
fprintf('Sum power constraint, decode 1 then 2:\nPower: (%.4f, %.4f)\nSum Rate: %.4f\n', [2-p2, p2, -fval]);

%% time sharing individual constraint, TDM
% lambda: time share parameter
[lambda, fval] = fminbnd(@(l) -weight2*l.*R(g22./l) - (1-l).*R(g11./(1-l)), 0, 1);
fprintf('Individual power constraint, TDM:\n%.4f of time user 1 transmits w/ power %.4f\n%.4f of time user 2 transmits w/ power %.4f\nSum Rate: %.4f\n\n',...
    [1-lambda, 1/(1-lambda), lambda, 1/lambda, -fval]);

%% time sharing individual constraint, NOMA
% grid search on lambda: fraction of time that Rxs decode 1 before 2
lambdas = 0.1:0.1:0.2;
p1s = zeros(1, length(lambdas));
p2s = zeros(1,length(lambdas));
sumrates = zeros(1,length(lambdas));
R12s = zeros(1,length(lambdas));
R21s = zeros(1,length(lambdas));

for idx = 1:length(lambdas)
    l = lambdas(idx);
    lb = 1-l;
    cvx_begin
        variable p(2) nonnegative % power vector for [1,2] case, the other case the powers are (1-l*p(1))/lb, (1-l*p(2))/lb
        variable rs(2) nonnegative % rate of the user to be decoded first r(1)=r121, r(2)=r212
        maximize l*rs(1) + l*weight2*log(1+g22*p(2)) + lb*log(1+g11*(1-l*p(1))/lb) + lb*weight2*rs(2)
        subject to
            p(1) <= 1/l; p(2) <= 1/l;
            (1+g12*p(2)).*(exp(rs(1))-1) <= g11*p(1);
    cvx_end
%     if g11*g22 > g12*g21 && g22>g12 
%         % when dec by [2,1], R2_at_1 is the bottle neck, rate sum R21
%         % independent of P1 -> so better set as zero and allocate all P1 to
%         % the other case, so p1 = [1/lambda, 0]
%         R21 = lb*R(g12*(1-l*p2)/lb);
%         R12 = 
%     end
end

    
% 
% 
% 
% 
% for idx = 1:length(lambdas)
%     l = lambdas(idx);
%     lb = 1-l;
%     cvx_begin quiet
%         variables p1 p2
%         maximize l*log(1+g22*p2 + g21*p1) + lb*log(1 + g11*(1-l*p1)/lb + g12*(1-l*p2)/lb)
%         subject to
%             p1 >= 0;
%             p2 >= 0;
%             p1 <= 1/l;
%             p2 <= 1/l;
%             r11 <= log(1+g11*p1);
%     cvx_end
%     p1s(idx) = p1;
%     p2s(idx) = p2;
%     sumrates(idx) = 0.5/log(2)*cvx_optval;
% end
%%
% get final result
[MaxRate, idx] = max(sumrates);
l = lambdas(idx);
fprintf('Individual power constraint, NOMA w/ time share:\n%.4f of time use power vector (%.4f, %.4f)\n%.4f of time use power vector (%.4f, %.4f)\nSum Rate: %.4f\n', ...
    [l, p1, p2, 1-l, (1-l*p1)/(1-l), (1-l*p2)/(1-l), MaxRate]);

%% individual power constraint (1,1), HK, no time sharing
% dec x1c -> x2c -> private
a2s = 0:0.001:1;
a1s = zeros(1, length(a2s));
sumrates = zeros(1, length(a2s));
p11 = 1;
p21 = 1;
R1 = @(a1,a2) R(g11*a1*p11./(g12*a2*p21 + 1)) ...
    + min(R(g21*(1-a1)*p11./(g21*a1*p11+g22*p21+1)), R(g11*(1-a1)*p11./(g11*a1*p11+g12*p21+1)));
R2 = @(a1,a2) R(g22*a2*p21./(g21*a1*p11 + 1)) ...
    + min(R(g12*(1-a2)*p21./(g11*a1*p11+g12*a2*p21+1)), R(g22*(1-a2)*p21./(g21*a1*p11+g22*a2*p21+1)));

weight2=1.1;
for idx = 1:length(a2s)
    a2 = a2s(idx);
    [a1s(idx), fval] = fminbnd(@(a) -R(g11*a*p11./(g12*a2*p21 + 1)) ...
        - weight2*min(R(g12*(1-a2)*p21./(g11*a*p11+g12*a2*p21+1)), R(g22*(1-a2)*p21./(g21*a*p11+g22*a2*p21+1))) ...
        - min(R(g21*(1-a)*p11./(g21*a*p11+g22*p21+1)), R(g11*(1-a)*p11./(g11*a*p11+g12*p21+1))) ...
        - weight2*R(g22*a2*p21./(g21*a*p11 + 1)), -1e-4, 1+1e-4);
    sumrates(idx) = -fval;
end
[MaxRate, idx] = max(sumrates);
fprintf('Individual power constraint, HK w/o time sharing, x1c->x2c->private:\n');
fprintf('Power on private message: (%.4f, %.4f)\n', [a1s(idx), a2s(idx)]);
fprintf('Individual rates: R1=%.4f, R2=%.4f\nSum Rate: %.4f\n\n', [R1(a1s(idx), a2s(idx)), R2(a1s(idx), a2s(idx)), MaxRate]);

% dec x2c -> x1c -> private
a2s = 0:0.001:1;
a1s = zeros(1, length(a2s));
sumrates = zeros(1, length(a2s));
p12 = 1;
p22 = 1;
R1 = @(a1,a2) R(g11*a1*p12./(g12*a2*p22 + 1)) ...
    + min(R(g21*(1-a1)*p12./(g21*a1*p12+g22*a2*p22+1)), R(g11*(1-a1)*p12./(g11*a1*p12+g12*a2*p22+1)));
R2 = @(a1,a2) R(g22*a2*p22./(g21*a1*p12 + 1)) ...
    + min(R(g12*(1-a2)*p22./(g11*p12+g12*a2*p22+1)), R(g22*(1-a2)*p22./(g21*p12+g22*a2*p22+1)));

for idx = 1:length(a2s)
    a2 = a2s(idx);
    [a1s(idx), fval] = fminbnd(@(a) -R(g11*a*p12./(g12*a2*p22 + 1)) ...
        - weight2*min(R(g12*(1-a2)*p22./(g11*p12+g12*a2*p22+1)), R(g22*(1-a2)*p22./(g21*p12+g22*a2*p22+1))) ...
        - min(R(g21*(1-a)*p12./(g21*a*p12+g22*a2*p22+1)), R(g11*(1-a)*p12./(g11*a*p12+g12*a2*p22+1))) ...
        - weight2*R(g22*a2*p22./(g21*a*p12 + 1)), -1e-4, 1+1e-4);
    sumrates(idx) = -fval;
end
[MaxRate, idx] = max(sumrates);
fprintf('Individual power constraint, HK, x2c->x1c->private:\n');
fprintf('Power on private message: (%.4f, %.4f)\n', [a1s(idx), a2s(idx)]);
fprintf('Individual rates: R1=%.4f, R2=%.4f\nSum Rate: %.4f\n\n', [R1(a1s(idx), a2s(idx)), R2(a1s(idx), a2s(idx)), MaxRate]);

%% individual power constraint (1,1), HK, with time sharing
% line search on time sharing factor
lambda = 0:0.1:1;

