%% non-time-share regions
BD_1_12_1 = @(E1, g11, g12, eb1) -1/g12 - g11/g12*E1/eb1;
BD_2_12_1 = @(E1, g21, g22, eb1) -1/g22 - g21/g22*E1/eb1;
BD_2_12_2 = @(g22, eb2) -eb2/g22;
BD_1_21_1 = @(g11, eb1) -eb1/g11;
BD_1_21_2 = @(E1, g11, g12, eb2) -eb2/g12*(1+g11*E1);
BD_2_21_2 = @(E1, g21, g22, eb2) -eb2/g22*(1+g21*E1);

H = [1, 0.8; 0.6, 1];
G = H.^2;
bvec = [0.4,0.3];
exp_b_vec = 1 - 2.^(2.*bvec);
Npt = 1000;
maxE = 10;
rangeE1 = linspace(0,maxE,Npt);
rangeE2 = linspace(0,maxE,Npt);
[E1, E2] = meshgrid(rangeE1, rangeE2);
% order [1,2], [1,2]
cond_1_12_1 = G(1,1)*E1 + exp_b_vec(1)*G(1,2)*E2 + exp_b_vec(1) >= 0;
cond_2_12_1 = G(2,1)*E1 + exp_b_vec(1)*G(2,2)*E2 + exp_b_vec(1) >= 0;
cond_2_12_2 = G(2,2)*E2 + exp_b_vec(2) >= 0;
cond_1_21_2 = exp_b_vec(2)*G(1,1)*E1 + G(1,2)*E2 + exp_b_vec(2) >= 0;
cond_1_21_1 = G(1,1)*E1 + exp_b_vec(1) >= 0;
cond_2_21_2 = exp_b_vec(2)*G(2,1)*E1 + G(2,2)*E2 + exp_b_vec(2) >= 0;
% dark red region
output_12_12 = zeros(Npt);
output_12_12(cond_1_12_1 & cond_2_12_1 & cond_2_12_2) = 1;
output_12_12 = cat(3, 0.5*output_12_12, zeros(Npt), zeros(Npt));
% dark green region
output_12_21 = zeros(Npt);
output_12_21(cond_1_12_1 & cond_2_21_2) = 1;
output_12_21 = cat(3, zeros(Npt), 0.5*output_12_21, zeros(Npt));
% dark blue region
output_21_12 = zeros(Npt);
output_21_12(cond_1_21_1 & cond_1_21_2 & cond_2_12_1 & cond_2_12_2) = 1;
output_21_12 = cat(3, zeros(Npt), zeros(Npt), output_21_12);
% dark yellow region
output_21_21 = zeros(Npt);
output_21_21(cond_1_21_2 & cond_1_21_1 & cond_2_21_2) =1;
output_21_21 = cat(3, 0.3*output_21_21, 0.3*output_21_21, zeros(Npt));

%imshow(1-output_12_12, 'xdata', rangeE1, 'ydata', rangeE2); hold on;
%imshow((output_12_12 + output_12_21 + output_21_12 + output_21_21), 'xdata', rangeE1, 'ydata', rangeE2);
background = min(~(output_12_12 + output_12_21 + output_21_21), [], 3);
background = cat(3, background, background, background);
imshow(background + (output_12_12 + output_12_21 + output_21_21), 'xdata', rangeE1, 'ydata', rangeE2);
set(gca,'YDir','normal'); hold on;
axis on;

%% TDM regions
t = linspace(0,1,500);
t = t(2:end-1);
E1_TDM = t.*(2.^(2*bvec(1)./t) - 1)/G(1,1);
E2_TDM = (1-t).*(2.^(2*bvec(2)./(1-t)) - 1)/G(2,2);
h(2) = plot(E1_TDM, E2_TDM, 'r-', 'LineWidth', 2);

% plot(rangeE1, BD_2_12_1(rangeE1, G(2,1), G(2,2), exp_b_vec(1)), 'r-', 'LineWidth', 2);
% plot(rangeE1, ones(1,Npt)*BD_2_12_2(G(2,2), exp_b_vec(2)), 'r--', 'LineWidth', 2)
% 
% plot(rangeE1, BD_1_12_1(rangeE1, G(1,1), G(1,2), exp_b_vec(1)), 'g-', 'LineWidth', 2);
% plot(rangeE1, BD_2_21_2(rangeE1, G(2,1), G(2,2), exp_b_vec(2)), 'g--', 'LineWidth', 2);
% 
% plot(rangeE1, BD_1_21_2(rangeE1, G(1,1), G(1,2), exp_b_vec(2)), 'y-', 'LineWidth', 2);
% plot(ones(1,Npt)*BD_1_21_1(G(1,1), exp_b_vec(1)), rangeE2, 'y--', 'LineWidth', 2);

pt_12_12 = [ -exp_b_vec(1)/G(2,1)*(1 + G(2,2)*BD_2_12_2(G(2,2), exp_b_vec(2))); BD_2_12_2(G(2,2), exp_b_vec(2))];
pt_12_21 = [G(1,1)/exp_b_vec(1), G(1,2); G(2,1), G(2,2)/exp_b_vec(2)]\[-1;-1];
pt_21_21 = [BD_1_21_1(G(1,1), exp_b_vec(1)); BD_1_21_2(BD_1_21_1(G(1,1), exp_b_vec(1)), G(1,1), G(1,2), exp_b_vec(2))];

rangeE1_12_12 = [pt_12_12(1), maxE];
rangeE2_21_21 = [pt_21_21(2), maxE];
h(1) = plot([pt_12_21(1), pt_21_21(1)], [pt_12_21(2), pt_21_21(2)], 'c-', 'LineWidth', 2);
plot([pt_12_21(1), pt_12_12(1)], [pt_12_21(2), pt_12_12(2)], 'c-', 'LineWidth', 2);
plot(rangeE1_12_12, [1,1]*pt_12_12(2), 'c-', 'LineWidth', 2);
plot([1,1]*pt_21_21(1), rangeE2_21_21, 'c-', 'LineWidth', 2);
hold off;
xlabel('E_1');
ylabel('E_2');
legend(h, {'Conv hull of different orderings', 'TDM'})

%% rigid time-share with 2 sub-strategies
H = [1, 0.2; 0.8, 1];
G = H.^2;
bvec = [0.1,0.1];
g11 = G(1,1); g12 = G(1,2); g21 = G(2,1); g22 = G(2,2);
% case 1+1 {[1,2], [1,2]} + {[1,2], [1,2]}
w = [1,1];
s2 = [1,1];
cvx_begin
cvx_solver mosek
    variable E1(2) nonnegative
    %variable E2(2) nonnegative
    variable s1(2) nonnegative
    %variable s2(2) nonnegative
    minimize w*E1
    subject to
        [g11*E1(1), 1+g12*E1(2); s1(1), 1] == semidefinite(2); % Rx 1 dec u1
        [g21*E1(1), 1+g22*E1(2); s1(1), 1] == semidefinite(2); % Rx 2 dec u2
        g22*E1(2) == s1(2);
        %[1+G(1,:)*E2, 1+g12*E2(2); s2(1), 1] == semidefinite(2);
        %[1+G(2,:)*E2, 1+g22*E2(2); s2(1), 1] == semidefinite(2);
        %1 + g22*E2(2) >= s2(2);
        log(1+s1(1)) >= 2*log(2)*bvec(1);
        log(1+s1(2)) >= 2*log(2)*bvec(2);
cvx_end
E1
0.5*(log2(s1))

% case 1+2
        
        