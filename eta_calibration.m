close all;
clear all;
clc;

%% MRR 4.6 data
L_ring = 178.98e-6;
R_ring = L_ring/(2*pi);
ng = 4.1850;
neff = 2.45;
L_b1 = 116.8e-6;
L_r1 = 47.12e-6;
L_b2 = 47.12e-6;
L_r2 = 47.12e-6;
k0 = 0.0441;
alpha = 800;
A = 1e9;

lambda0 = 1550.391e-9;
f0 = mrr_asym.c / lambda0;
w0 = 2*pi*f0;

loss_dB_rt = alpha * L_ring;
transmission_rt = 10^(-loss_dB_rt/10);
etha = 1 - transmission_rt;

ref = mrr_asym(R_ring, neff, ng, k0, L_b1, L_r1, L_b2, L_r2, ng, ng, A, alpha);
Qi = ref.Q_Wu(f0, ng, etha);

% set 1
P1_data = [0, 1.52, 2.93];
a0_set1 = [8.0127e10, 7.6654e10, 7.0308e10];
b0_set1 = [8.1233e9,  1.1596e10, 1.7942e10];

% set 2
P2_data = [0, 2.92, 5.75];
a0_set2 = [8.0127e10, 7.8875e10, 7.5318e10];
b0_set2 = [8.1233e9,  6.8712e9,  3.3143e9];

dn1 = zeros(size(P1_data));
for i = 1:numel(P1_data)
    Qe1_target = w0 / (a0_set1(i) - b0_set1(i));
    k1_target  = 1 - exp(-w0*ng*L_ring/(mrr_asym.c*Qe1_target));
    nb1_sol = invert_kappa(ref, lambda0, 1, k1_target, ng, ng-0.02, ng+0.02);
    dn1(i) = nb1_sol - ng;
end

[~, k2_ng] = ref.kappa(lambda0); 
Qe2_ng = ref.Q_Wu(f0, ng, k2_ng);
S_baseline = a0_set2(1) + b0_set2(1);
Qi_cal = 1 / (S_baseline/w0 - 1/Qe2_ng);

dn2 = zeros(size(P2_data));
dn2(1) = 0;  % P2=0 -> nb2=ng
for i = 2:numel(P2_data)
    Qe2_target = 1 / ((a0_set2(i)+b0_set2(i))/w0 - 1/Qi_cal);
    k2_target  = 1 - exp(-w0*ng*L_ring/(mrr_asym.c*Qe2_target));
    % dn/dT>0 for Si -> always adopt positive solution.
    nb2_sol = invert_kappa(ref, lambda0, 2, k2_target, ng, ng, ng+0.02);
    dn2(i) = nb2_sol - ng;
end

%% linear fit (tuning efficiency)
p1fit = polyfit(P1_data, dn1, 1);
p2fit = polyfit(P2_data, dn2, 1);
 
eta1 = p1fit(1);   % 1/mW
eta2 = p2fit(1);   % 1/mW
 
fprintf('eta1 = %.4e 1/mW  (intercept %.4e)\n', eta1, p1fit(2));
fprintf('eta2 = %.4e 1/mW  (intercept %.4e)\n', eta2, p2fit(2));
 
%% Visualization
figure('Color', 'w', 'Position', [100 100 1000 450]);
 
subplot(1,2,1);
plot(P1_data, dn1, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor','k', 'MarkerSize', 8); hold on;
Pfit1 = linspace(0, max(P1_data), 100);
plot(Pfit1, polyval(p1fit, Pfit1), 'k--', 'LineWidth', 1.5);
xlabel('P_1  [mW]'); ylabel('\Delta n_{b1}');
title(sprintf('Heater 1:  \\eta_1 = %.3e  /mW', eta1));
legend('Data Test 1', 'Linear fit', 'Location', 'northwest');
grid on;
 
subplot(1,2,2);
plot(P2_data, dn2, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor','k', 'MarkerSize', 8); hold on;
Pfit2 = linspace(0, max(P2_data), 100);
plot(Pfit2, polyval(p2fit, Pfit2), 'k--', 'LineWidth', 1.5);
xlabel('P_2  [mW]'); ylabel('\Delta n_{b2}');
title(sprintf('Heater 2:  \\eta_2 = %.3e  /mW', eta2));
legend('Data Test 2', 'Linear fit', 'Location', 'northwest');
grid on;
 
sgtitle('Tuning Efficiency of TiN microheaters');
 
 
function nb_sol = invert_kappa(ref_mrr, lambda0, which, k_target, ng, nb_lo, nb_hi)
    grid = linspace(nb_lo, nb_hi, 4001);
    f = zeros(size(grid));
    for i = 1:numel(grid)
        f(i) = local_k(ref_mrr, lambda0, which, grid(i), ng) - k_target;
    end
 
    idx = find(f(1:end-1).*f(2:end) < 0);
 
    roots = zeros(size(idx));
    for j = 1:numel(idx)
        fun = @(nb) local_k(ref_mrr, lambda0, which, nb, ng) - k_target;
        roots(j) = fzero(fun, [grid(idx(j)), grid(idx(j)+1)]);
    end
 
    [~, best] = min(abs(roots - ng));
    nb_sol = roots(best);
end
 
function k = local_k(ref_mrr, lambda0, which, nb, ng)
    trial = ref_mrr;
    if which == 1
        trial.nb1 = nb; trial.nb2 = ng;
        [k, ~] = trial.kappa(lambda0);
    else
        trial.nb1 = ng; trial.nb2 = nb;
        [~, k] = trial.kappa(lambda0);
    end
end








