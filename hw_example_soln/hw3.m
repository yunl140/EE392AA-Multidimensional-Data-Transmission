clc; close all;
%% 4.13
clear;
disp('------------ 4.13 ------------')
b = [2 4 6];
dfree = [10 6 5 1];
r = [1/2 2/3 3/4 1];
Ntones = 48;
snr = 10^1.45;
symbol_rate = .25;
% part a
disp('------ part a ------')
possible_rates = symbol_rate*Ntones*b'*r
% part b
disp('------ part b ------')
snr = 10^(1.45+.63);
feasibleMCS = log10(2*qfunc(sqrt(3*dfree'*snr*(2.^b-1).^(-1))))'<-6.9;
[maxrate, idx]=max(possible_rates.*feasibleMCS, [], 'all', 'linear');
[M_idx, r_idx]=ind2sub([3 4],idx);
fprintf('max rate: %d Mbps, %.2f bits/dimension, achieved at %d QAM, code rate %.2f\n\n', maxrate, maxrate/Ntones/symbol_rate/2, 2^b(M_idx), r(r_idx));
% part c
disp('------ part c ------')
snr = 10^(1.45+.63-.3);
feasibleMCS = log10(2*qfunc(sqrt(3*dfree'*snr*(2.^b-1).^(-1))))'<-6.9;
[maxrate, idx]=max(possible_rates.*feasibleMCS, [], 'all', 'linear');
[M_idx, r_idx]=ind2sub([3 4],idx);
fprintf('max rate: %d Mbps, %.2f bits/dimension, achieved at %d QAM, code rate %.2f\n\n', maxrate, maxrate/Ntones/symbol_rate/2, 2^b(M_idx), r(r_idx));
% part d
disp('------ part d ------')
snr = 10^(1.45+.63+.55);
b = [2 4 6 8];
possible_rates = symbol_rate*Ntones*b'*r;
possible_bit_per_dim = b'*r/2;
feasibleMCS = log10(2*qfunc(sqrt(3*dfree'*snr*(2.^b-1).^(-1))))'<-6.9;
[maxrate, idx]=max(possible_rates.*feasibleMCS, [], 'all', 'linear');
[M_idx, r_idx]=ind2sub([4 4],idx);
fprintf('max rate: %d Mbps, %.2f bits/dimension, achieved at %d QAM, code rate %.2f\n\n', maxrate, maxrate/Ntones/symbol_rate/2, 2^b(M_idx), r(r_idx));
% part e
disp('------ part e ------')
ratio = (3/4*6)/(1/2*1)
fprintf('\n');

%% P4.14
clear;
disp('------------ 4.14 ------------')
gbounds = [0 .0105 .0223 .0357 .0511 .0693 .0916 .1204 .1609 .2303 .4];
gs = (gbounds(2:end)+gbounds(1:end-1))/2;
EsN0 = 10^4.3;
%EsN0 = 10^3.3;
% part a
disp('------ part a ------')
ave_g = mean(gs);
ave_snr = EsN0*ave_g;
ave_snr_db = 10*log10(ave_snr);
fprintf('<g>: %.4f, <SNR>: %.4f dB\n', ave_g, ave_snr_db);
M = [4 16 64];
% if consider the average snr
Pebar = 2*qfunc(sqrt(3*1*ave_snr./(M-1)));
% if consider the worst channel
Pebar = 2*qfunc(sqrt(3*1*gs(1)*EsN0./(M-1)));

% part b
disp('------ part b ------')
gaps_db = [9 3];
gaps = 10.^([.9 .3]);
for gap_idx = 1:2
    gap = gaps(gap_idx);
    for m = M
        for numG = 10:-1:1
            log_Kma_over_gap = 1/numG*(10*log2(m)-sum(log2(gs(11-numG:end))));
            if -log2(gs(11-numG))<log_Kma_over_gap
                Es = gap*(2^log_Kma_over_gap-1./gs(11-numG:end));
                meanEs = sum(Es)/10;
                margin=10*log10(EsN0/meanEs);
                fprintf('Gap: %d dB, M: %d, |Gstar|: %d, margin: %.4f dB\n', gaps_db(gap_idx), m, numG, margin);
                break;
            end
        end
    end
end

% part c
disp('------ part c ------')
for gap_idx = 1:2
    gap = gaps(gap_idx);
    for numG = 10:-1:1
        Kra = (10*EsN0+gap*sum(1./gs(11-numG:end)))/numG;
        if Kra>gap/gs(11-numG)
            Es = Kra - gap./gs(11-numG:end);
            meanEs = sum(Es)/10;
            b = 0.1*sum(log2(Kra*gs(11-numG:end)/gap));
            fprintf('Gap: %d dB, |Gstar|: %d, Echeck: %.4f, b: %.4f\n\n', gaps_db(gap_idx), numG, meanEs, b);
            break;
        end
    end
end

%% 4.22
clear;
disp('------------ 4.22 ------------')
% part a
disp('------ part a ------')
gn = [20 30 40];
gap=10^.9;
Ebase = 2*gap./gn*3;
incremental_table = 4.^(0:2)'*Ebase
extend_table = [incremental_table,incremental_table(:,end)];
% part b
disp('------ part b ------')
datarate = (4*2+2*2)*10;
fprintf('Data rate is %d Mbps\n\n', datarate);
% part c
disp('------ part c ------')
Euse = sum(extend_table(1,:))+extend_table(2,end);
margin=10*log10(16/Euse);
fprintf('max margin is %.2f dB\n\n', margin);

%% 4.16
clear;
disp('------------ 4.16 ------------')
% part a
g_func = @(M,dinv) qfuncinv(1e-5)^2/3*(M-1)'*dinv;
% part b
disp('------ part b ------')
r=[1/2 2/3 3/4 5/6 7/8];
dfree = [10 6 5 4 3];
symbol_rate = .25;
Ntones = 48;
n = [1 2 3 4];
target_gs = g_func(4.^n, dfree.^(-1));
target_gs_snr = 10*log10(target_gs)
% part c
disp('------ part c ------')
data_rates = Ntones*symbol_rate*2*n'*r
% part d
disp('------ part d ------')
target_gs_lin = [0 reshape(target_gs', 1, []) Inf];
pg = exp(-target_gs_lin(1:end-1)/10) - exp(-target_gs_lin(2:end)/10)
% part e
disp('------ part e ------')
target_gs_lin(end) = target_gs_lin(end-1);
g_center = (target_gs_lin(1:end-1)+target_gs_lin(2:end))/2;
ave_rate = pg*log2(1+g_center)';
fprintf('using middle point: %.4f bits/dim, %.4f Mbps\n', ave_rate, Ntones*symbol_rate*ave_rate);
ave_rate = pg*log2(1+target_gs_lin(1:end-1))';
fprintf('using min point: %.4f bits/dim, %.4f Mbps\n\n', ave_rate, Ntones*symbol_rate*ave_rate);
% part f
N = 48;
g_sample = -10*log(1-rand(1,N));
max_sample=max(g_sample);
stop_bin = find(target_gs_lin>max_sample,1);
figure(1); hold on;
histogram(g_sample, target_gs_lin(1:stop_bin), 'Normalization','probability');
plot((target_gs_lin(1:stop_bin)+target_gs_lin(2:stop_bin+1))/2, pg(1:stop_bin));
hold off;
% part g
N = 20;
g_sample = -10*log(1-rand(1,N));
max_sample=max(g_sample);
stop_bin = find(target_gs_lin>max_sample,1);
edges = target_gs_lin(1:stop_bin);
figure(2); hold on;
histogram(g_sample, target_gs_lin(1:stop_bin), 'Normalization','probability');
plot((target_gs_lin(1:stop_bin)+target_gs_lin(2:stop_bin+1))/2, pg(1:stop_bin));
hold off;

%% 4.15
clear;
disp('------------ 4.15 ------------')
% part a
disp('------ part a ------')
N0 = -174+6;
band_db = 10*log10([20e6, 40e6]);
fprintf('noise power is %.2f dBm (20MHz); %.2f dBm (40MHz)\n', ...
    N0+band_db(1), N0+band_db(2));
fprintf('the PSD is %d dBm/Hz for both options\n\n', N0);
% part b
disp('------ part b ------')
snr_db = 15-band_db-90-N0;
fprintf('SNR is %.2f dB (20MHz); %.2f dB (40MHz)\n\n', snr_db(1), snr_db(2));
% part c
disp('------ part c ------')
Ntones = [108 48];
fprintf('ratio is %.2f\n\n', Ntones(1)/Ntones(2));
% part d
disp('------ part d ------')
pg = [.25 .35 .20 .15 .05];
gs = [0 3 6 12 24];
ave_g = -90-N0+10*log10(2);
gtilde = ave_g-10*log10(pg*gs');
fprintf('gtilde is %.2f dB\n\n', gtilde);
% part e
disp('------ part e ------')
Pe_func = @(M, d, g) 2*qfunc(sqrt(3*d'*g/(M-1)));
r=[1/2 2/3];
dfree = [10 6];
ns = [1 2 3 4];

% the first bin is unusable -> Pout >=0.25 -> r <=3/4
% if the second bin is also off -> Pout >= 0.6, can't be recovered by any
% code
% so, all 4 nonzero bins are on in RA loading, r in [1/2 2/3 3/4]
pg = pg(2:end);
% original
%gtilde_norm = 10.^((15-band_db-10*log10(2)+gtilde)/10)/.75; % Normalize to Ex_bar=1
% no .75 increase
gtilde_norm = 10.^((15-band_db-10*log10(2)+gtilde)/10); % Normalize to Ex_bar=1
% match John's
%gtilde=81-10*log10(5.3);
%gtilde_norm = 10.^(([-60 -63]+gtilde)/10); % Normalize to Ex_bar=1
max_rate = zeros(1,2);
opt_code_rate = zeros(1,2);
opt_n = zeros(1,2);

for band_idx = 1:2
    gs = [3 6 12 24]*gtilde_norm(band_idx);
    for n = ns
        Pe = Pe_func(4^n, dfree, gs);
        %feasible_r = Pe*pg' < 1e-6;
        feasible_r = Pe(:,1) < 1e-6;
        if any(feasible_r)
            r_tmp = r(find(feasible_r, 1, 'last'));
            rate = 2*r_tmp*n;
            if max_rate(band_idx) < rate
               max_rate(band_idx) = rate;
               opt_code_rate(band_idx) = r_tmp;
               opt_n(band_idx) = n;
            end
        end
    end
end
fprintf('20MHz: %.2f bits/tone, r=%.2f, %d QAM\n40MHz: %.2f bits/tone, r=%.2f, %d QAM\n\n',...
    max_rate(1), opt_code_rate(1), 4^opt_n(1), max_rate(2), opt_code_rate(2), 4^opt_n(2));
% part f
disp('------ part f ------')
fprintf('data rate: %d Mbps (20MHz), %d Mbps (40MHz)\n', max_rate(1)*48/4, max_rate(2)*108/4);