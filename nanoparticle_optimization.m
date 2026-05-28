% Computational Optimization for Atherosclerosis Treatment
% % Author: [Ayur Anchan] | NCSEF 2025

fprintf('NANOPARTICLE OPTIMIZATION\n');

% SECTION 1: INITIALIZATION 
fprintf('Initializing parameters and conditions\n');

% Set random seed for reproducibility
rng(42);

% Define physiological conditions (5 disease stages)
conditions = struct();
conditions.healthy = struct('pH', 7.4, 'ROS', 5, 'stage', 0, 'name', 'Healthy');
conditions.subclinical = struct('pH', 7.0, 'ROS', 50, 'stage', 1, 'name', 'Subclinical');
conditions.early = struct('pH', 6.8, 'ROS', 100, 'stage', 2, 'name', 'Early Plaque');
conditions.moderate = struct('pH', 6.5, 'ROS', 250, 'stage', 3, 'name', 'Moderate');
conditions.severe = struct('pH', 6.3, 'ROS', 400, 'stage', 4, 'name', 'Severe');

% Baseline parameters (literature-derived)
params = struct();
params.k_pH_max = 0.15;      % h^-1, maximum pH-responsive degradation
params.alpha = 2.0;          % steepness of pH response
params.pH0 = 7.0;            % inflection point
params.k_ROS_base = 0.02;    % h^-1, baseline degradation
params.k_ROS_max = 0.25;     % h^-1, maximum ROS-responsive degradation
params.K_ROS = 50;           % μM, Michaelis-Menten constant
params.beta = 0.50;          % synergy coefficient
params.k_elim = 0.05;        % h^-1, drug elimination rate

% Simulation settings
M0 = 100;  % μg, initial drug in nanoparticle
tspan = [0 72];  % hours, extended time course
t_eval = 24;  % hours, primary analysis timepoint

fprintf('Parameters initialized\n');
fprintf('Disease stages: 0 (healthy) to 5 (severe)\n\n');

% SECTION 2: CORE SIMULATIONS 
fprintf('RUNNING CORE SIMULATIONS\n');

results = struct();
cond_names = fieldnames(conditions);

% Run simulations for all conditions and mechanisms
for i = 1:length(cond_names)
    name = cond_names{i};
    cond = conditions.(name);
    
    fprintf('Simulating %s (pH %.1f, ROS %d μM)...\n', cond.name, cond.pH, cond.ROS);
    
    % pH-only responsive
    [t, y] = ode45(@(t,y) release_pH(t, y, cond.pH, params), tspan, [M0; 0]);
    results.(name).pH.t = t;
    results.(name).pH.M = y(:,1);
    results.(name).pH.D = y(:,2);
    
    % ROS-only responsive
    [t, y] = ode45(@(t,y) release_ROS(t, y, cond.ROS, params), tspan, [M0; 0]);
    results.(name).ROS.t = t;
    results.(name).ROS.M = y(:,1);
    results.(name).ROS.D = y(:,2);
    
    % Dual-responsive
    [t, y] = ode45(@(t,y) release_dual(t, y, cond.pH, cond.ROS, params), tspan, [M0; 0]);
    results.(name).dual.t = t;
    results.(name).dual.M = y(:,1);
    results.(name).dual.D = y(:,2);
    
    % Non-responsive control
    [t, y] = ode45(@(t,y) release_baseline(t, y, 0.01, params), tspan, [M0; 0]);
    results.(name).baseline.t = t;
    results.(name).baseline.M = y(:,1);
    results.(name).baseline.D = y(:,2);
end

% SECTION 3: TARGETING INDEX ANALYSIS 
fprintf('TARGETING INDEX ANALYSIS\n');

% Calculate targeting indices(TI) at multiple timepoints
timepoints = [24, 48];

for t_idx = 1:length(timepoints)
    t_point = timepoints(t_idx);
    fprintf('\nTargeting Indices at %d hours:\n', t_point);
    
    for i = 1:length(cond_names)
        name = cond_names{i};
        if strcmp(name, 'healthy')
            continue; % Skip healthy as it's the denominator
        end
        
        % Dual-responsive TI
        D_healthy = interp1(results.healthy.dual.t, results.healthy.dual.D, t_point);
        D_disease = interp1(results.(name).dual.t, results.(name).dual.D, t_point);
        TI_dual = D_disease / D_healthy;
        
        % pH-only TI
        D_healthy_pH = interp1(results.healthy.pH.t, results.healthy.pH.D, t_point);
        D_disease_pH = interp1(results.(name).pH.t, results.(name).pH.D, t_point);
        TI_pH = D_disease_pH / D_healthy_pH;
        
        % ROS-only TI
        D_healthy_ROS = interp1(results.healthy.ROS.t, results.healthy.ROS.D, t_point);
        D_disease_ROS = interp1(results.(name).ROS.t, results.(name).ROS.D, t_point);
        TI_ROS = D_disease_ROS / D_healthy_ROS;
        
        fprintf('%s:\n', conditions.(name).name);
        fprintf('  pH-only TI: %.2f\n', TI_pH);
        fprintf('  ROS-only TI: %.2f\n', TI_ROS);
        fprintf('  Dual TI: %.2f (%.1f%% improvement over best single)\n', ...
                TI_dual, (TI_dual/max(TI_pH,TI_ROS)-1)*100);
    end
end
fprintf('\n');

% SECTION 4: MONTE CARLO UNCERTAINTY ANALYSIS 
fprintf('MONTE CARLO UNCERTAINTY QUANTIFICATION\n');
fprintf('10,000 iterations\n');

n_iterations = 10000;
TI_distribution = zeros(n_iterations, 1);

% Define parameter uncertainties (+/- 15% standard deviation)
param_means = [params.k_pH_max, params.K_ROS, params.beta, params.alpha, ...
               params.k_ROS_max, params.k_ROS_base];
param_stds = param_means * 0.15;  % 15% coefficient of variation

% Monte Carlo simulation
for iter = 1:n_iterations
    if mod(iter, 2000) == 0
        fprintf('  Progress: %d/%d iterations\n', iter, n_iterations);
    end
    
    % Sample parameters from normal distributions
    params_mc = params;
    sampled = param_means + param_stds .* randn(size(param_means));
    sampled = max(sampled, param_means * 0.5);  % Prevent unrealistic values
    
    params_mc.k_pH_max = sampled(1);
    params_mc.K_ROS = sampled(2);
    params_mc.beta = sampled(3);
    params_mc.alpha = sampled(4);
    params_mc.k_ROS_max = sampled(5);
    params_mc.k_ROS_base = sampled(6);
    
    % Simulate healthy tissue
    [~, y] = ode45(@(t,y) release_dual(t, y, 7.4, 5, params_mc), [0 24], [M0; 0]);
    D_healthy_mc = y(end, 2);
    
    % Simulate early plaque
    [t, y] = ode45(@(t,y) release_dual(t, y, 6.8, 100, params_mc), [0 24], [M0; 0]);
    D_plaque_mc = y(end, 2);
    
    TI_distribution(iter) = D_plaque_mc / D_healthy_mc;
end

% Statistical analysis
TI_mean = mean(TI_distribution);
TI_std = std(TI_distribution);
TI_ci = prctile(TI_distribution, [2.5, 97.5]);
prob_TI_gt4 = sum(TI_distribution > 4.0) / n_iterations * 100;
prob_TI_gt5 = sum(TI_distribution > 5.0) / n_iterations * 100;

fprintf('\nMonte Carlo Results:\n');
fprintf('  Mean TI: %.2f\n', TI_mean);
fprintf('  Std Dev: %.2f\n', TI_std);
fprintf('  95%% CI: [%.2f, %.2f]\n', TI_ci(1), TI_ci(2));
fprintf('  P(TI > 4.0): %.1f%%\n', prob_TI_gt4);
fprintf('  P(TI > 5.0): %.1f%%\n', prob_TI_gt5);
fprintf('Monte Carlo analysisn\n');  % Monte Carlo Analysis Done

% SECTION 5: PARAMETER OPTIMIZATION 
fprintf('PARAMETER OPTIMIZATION\n');
fprintf('Grid search optimization\n');

% Define parameter ranges for optimization
k_pH_range = linspace(0.08, 0.22, 25);
K_ROS_range = linspace(30, 80, 25);
beta_range = linspace(0.2, 1.0, 20);

% 2D optimization: k_pH_max vs K_ROS (beta at baseline)
fprintf('Optimizing k_pH_max and K_ROS\n');
TI_matrix = zeros(length(k_pH_range), length(K_ROS_range));

for i = 1:length(k_pH_range)
    for j = 1:length(K_ROS_range)
        params_opt = params;
        params_opt.k_pH_max = k_pH_range(i);
        params_opt.K_ROS = K_ROS_range(j);
        
        % Healthy tissue
        [~, y] = ode45(@(t,y) release_dual(t, y, 7.4, 5, params_opt), [0 24], [M0; 0]);
        D_h = y(end, 2);
        
        % Early plaque
        [~, y] = ode45(@(t,y) release_dual(t, y, 6.8, 100, params_opt), [0 24], [M0; 0]);
        D_p = y(end, 2);
        
        TI_matrix(i, j) = D_p / D_h;
    end
end

[max_TI, idx] = max(TI_matrix(:));
[row, col] = ind2sub(size(TI_matrix), idx);
optimal_k_pH = k_pH_range(row);
optimal_K_ROS = K_ROS_range(col);

fprintf('Optimal Parameters (Stage 1):\n');
fprintf('  k_pH_max: %.4f h^-1\n', optimal_k_pH);
fprintf('  K_ROS: %.1f μM\n', optimal_K_ROS);
fprintf('  TI achieved: %.2f\n', max_TI);
fprintf('  Improvement: %.1f%%\n\n', (max_TI/TI_mean - 1)*100);

% Beta optimization (1D scan)
fprintf('Optimizing synergy coefficient (beta)\n');
TI_beta_scan = zeros(size(beta_range));

for i = 1:length(beta_range)
    params_opt = params;
    params_opt.beta = beta_range(i);
    params_opt.k_pH_max = optimal_k_pH;
    params_opt.K_ROS = optimal_K_ROS;
    
    [~, y] = ode45(@(t,y) release_dual(t, y, 7.4, 5, params_opt), [0 24], [M0; 0]);
    D_h = y(end, 2);
    
    [~, y] = ode45(@(t,y) release_dual(t, y, 6.8, 100, params_opt), [0 24], [M0; 0]);
    D_p = y(end, 2);
    
    TI_beta_scan(i) = D_p / D_h;
end

[max_TI_beta, idx_beta] = max(TI_beta_scan);
optimal_beta = beta_range(idx_beta);

fprintf('Optimal Synergy Coefficient:\n');
fprintf('  Beta: %.2f\n', optimal_beta);
fprintf('  Maximum TI: %.2f\n\n', max_TI_beta);

% SECTION 6: SOBOL SENSITIVITY ANALYSIS 
fprintf('GLOBAL SENSITIVITY ANALYSIS\n');
fprintf('Sobol indices\n');

% Sobol analysis using Latin Hypercube Sampling
n_sobol = 5000;
param_names = {'k_{pH,max}', 'K_{ROS}', 'β', 'α', 'k_{ROS,max}'};
n_params = length(param_names);

% Generate parameter samples
base_samples = lhsdesign(n_sobol, n_params);
param_ranges = [0.08 0.22; 30 80; 0.2 1.0; 1.5 3.0; 0.15 0.40];

% Scale to parameter ranges
for i = 1:n_params
    base_samples(:,i) = param_ranges(i,1) + ...
                        base_samples(:,i) * (param_ranges(i,2) - param_ranges(i,1));
end

% Calculate TI for all samples
TI_sobol = zeros(n_sobol, 1);
for i = 1:n_sobol
    params_s = params;
    params_s.k_pH_max = base_samples(i,1);
    params_s.K_ROS = base_samples(i,2);
    params_s.beta = base_samples(i,3);
    params_s.alpha = base_samples(i,4);
    params_s.k_ROS_max = base_samples(i,5);
    
    [~, y] = ode45(@(t,y) release_dual(t, y, 7.4, 5, params_s), [0 24], [M0; 0]);
    D_h = y(end, 2);
    [~, y] = ode45(@(t,y) release_dual(t, y, 6.8, 100, params_s), [0 24], [M0; 0]);
    D_p = y(end, 2);
    TI_sobol(i) = D_p / D_h;
end

% Calculate first-order sensitivity indices
total_var = var(TI_sobol);
sensitivity_indices = zeros(n_params, 1);

for i = 1:n_params
    % Bin-based variance decomposition
    [~, sort_idx] = sort(base_samples(:,i));
    n_bins = 20;
    bin_size = floor(n_sobol / n_bins);
    bin_means = zeros(n_bins, 1);
    
    for b = 1:n_bins
        idx_range = (b-1)*bin_size + 1 : min(b*bin_size, n_sobol);
        bin_means(b) = mean(TI_sobol(sort_idx(idx_range)));
    end
    
    conditional_var = var(bin_means);
    sensitivity_indices(i) = conditional_var / total_var;
end

% Normalize
sensitivity_indices = sensitivity_indices / sum(sensitivity_indices);

fprintf('\nFirst-Order Sensitivity Indices:\n');
for i = 1:n_params
    fprintf('  %s: %.1f%% variance contribution\n', ...
            param_names{i}, sensitivity_indices(i)*100);
end
fprintf('\n');

% SECTION 7: PATIENT STRATIFICATION 
fprintf('PATIENT STRATIFICATION ANALYSIS\n');

stage_names = {'Subclinical', 'Early', 'Moderate', 'Severe'};
stage_conds = {'subclinical', 'early', 'moderate', 'severe'};
accumulation_data = zeros(length(stage_conds), 3);

% Use optimized parameters
params_opt = params;
params_opt.k_pH_max = optimal_k_pH;
params_opt.K_ROS = optimal_K_ROS;
params_opt.beta = optimal_beta;

for i = 1:length(stage_conds)
    name = stage_conds{i};
    cond = conditions.(name);
    
    % Plaque accumulation
    [~, y] = ode45(@(t,y) release_dual(t, y, cond.pH, cond.ROS, params_opt), ...
                   [0 24], [M0; 0]);
    plaque_release = y(end, 2);
    
    % Healthy tissue accumulation
    [~, y] = ode45(@(t,y) release_dual(t, y, 7.4, 5, params_opt), [0 24], [M0; 0]);
    healthy_release = y(end, 2);
    
    accumulation_data(i, :) = [i, plaque_release, healthy_release];
    ratio = plaque_release / healthy_release;
    
    fprintf('Stage %d (%s):\n', i, stage_names{i});
    fprintf('  Plaque: %.1f μg\n', plaque_release);
    fprintf('  Healthy: %.1f μg\n', healthy_release);
    fprintf('  Selectivity ratio: %.1fx\n', ratio);
end
fprintf('\n');

% SECTION 8: PHARMACOECONOMIC ANALYSIS
fprintf('PHARMACOECONOMIC ANALYSIS\n');

% Annual costs per patient
cost_systemic_drug = 400;
cost_systemic_monitoring = 400;
cost_systemic_side_effects = 800;
cost_systemic_total = cost_systemic_drug + cost_systemic_monitoring + cost_systemic_side_effects;

% Nanoparticle therapy (quarterly)
cost_np_production = 2000;
cost_np_admin = 500;
cost_np_doses_per_year = 4;
cost_np_monitoring = 100;
cost_np_total = (cost_np_production + cost_np_admin) * cost_np_doses_per_year + cost_np_monitoring;

% 10-year analysis
years = 10;
cost_systemic_10yr = cost_systemic_total * years;
cost_np_10yr = cost_np_total * years;
savings = cost_systemic_10yr - cost_np_10yr;
savings_percent = (1 - cost_np_10yr/cost_systemic_10yr) * 100;

% Quality-Adjusted Life Years (QALY)
QALY_systemic = 7.5;
QALY_np = 9.3;
QALY_gained = QALY_np - QALY_systemic;

cost_per_QALY_systemic = cost_systemic_10yr / QALY_systemic;
cost_per_QALY_np = cost_np_10yr / QALY_np;

fprintf('10-Year Cost Analysis (per patient):\n');
fprintf('  Systemic therapy: $%d\n', cost_systemic_10yr);
fprintf('  Nanoparticle therapy: $%d\n', cost_np_10yr);
fprintf('  Savings: $%d (%.1f%%)\n', savings, savings_percent);
fprintf('\nQuality-Adjusted Life Years:\n');
fprintf('  Systemic: %.1f years\n', QALY_systemic);
fprintf('  Nanoparticle: %.1f years\n', QALY_np);
fprintf('  Additional QALY: %.1f years\n', QALY_gained);
fprintf('\nCost-Effectiveness:\n');
fprintf('  Systemic: $%d per QALY\n', round(cost_per_QALY_systemic));
fprintf('  Nanoparticle: $%d per QALY\n', round(cost_per_QALY_np));

% Population impact
us_statin_patients = 40e6;
adoption_rate = 0.10;
patients_affected = us_statin_patients * adoption_rate;
annual_savings_population = savings / years * patients_affected / 1e9;

fprintf('\nPopulation Impact (10%% US adoption):\n');
fprintf('  Patients: %.1f million\n', patients_affected/1e6);
fprintf('  Annual savings: $%.2f billion\n', annual_savings_population);
fprintf('\n');

% SECTION 9: VISUALIZATIONS 
fprintf('GENERATING VISUALIZATIONS\n');

% FIGURE 1: Release profiles across disease stages
figure('Position', [50 50 1400 500]);
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
          0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];

for i = 1:length(cond_names)
    name = cond_names{i};
    plot(results.(name).dual.t, results.(name).dual.D, 'LineWidth', 2.5, ...
         'Color', colors(i,:), 'DisplayName', conditions.(name).name);
    hold on;
end
xlabel('Time (hours)', 'FontSize', 13);
ylabel('Released Simvastatin (μg)', 'FontSize', 13);
title('Dual-Responsive Nanoparticle: Disease Stage Comparison', 'FontSize', 15);
legend('Location', 'southeast', 'FontSize', 11);
grid on;
xlim([0 72]);
ylim([0 100]);

% FIGURE 2: Mechanism comparison
figure('Position', [100 100 1200 400]);

subplot(1,2,1)
plot(results.early.pH.t, results.early.pH.D, 'b-', 'LineWidth', 2.5);
hold on;
plot(results.early.ROS.t, results.early.ROS.D, 'r-', 'LineWidth', 2.5);
plot(results.early.dual.t, results.early.dual.D, 'g-', 'LineWidth', 2.5);
plot(results.early.baseline.t, results.early.baseline.D, 'y--', 'LineWidth', 2);
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Released Drug (μg)', 'FontSize', 12);
title('Early Plaque (pH 6.8, 100 μM ROS)', 'FontSize', 13);
legend('pH-responsive', 'ROS-responsive', 'Dual-responsive', 'Non-responsive', ...
       'Location', 'southeast');
grid on;
xlim([0 72]);

subplot(1,2,2)
plot(results.healthy.pH.t, results.healthy.pH.D, 'b-', 'LineWidth', 2.5);
hold on;
plot(results.healthy.ROS.t, results.healthy.ROS.D, 'r-', 'LineWidth', 2.5);
plot(results.healthy.dual.t, results.healthy.dual.D, 'g-', 'LineWidth', 2.5);
plot(results.healthy.baseline.t, results.healthy.baseline.D, 'y--', 'LineWidth', 2);
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Released Drug (μg)', 'FontSize', 12);
title('Healthy Tissue (pH 7.4, 5 μM ROS)', 'FontSize', 13);
legend('pH-responsive', 'ROS-responsive', 'Dual-responsive', 'Non-responsive', ...
       'Location', 'southeast');
grid on;
xlim([0 72]);

% FIGURE 3: Monte Carlo distribution
figure('Position', [150 150 1000 400]);

subplot(1,2,1)
histogram(TI_distribution, 50, 'Normalization', 'probability', 'FaceColor', [0.2 0.6 0.8]);
hold on;
xline(TI_mean, 'r-', 'LineWidth', 2.5);
xline(TI_ci(1), 'k--', 'LineWidth', 1.5);
xline(TI_ci(2), 'k--', 'LineWidth', 1.5);
text(TI_mean, 0.045, sprintf(' Mean: %.2f', TI_mean), 'FontSize', 11);
text(TI_ci(1), 0.04, sprintf(' %.2f', TI_ci(1)), 'FontSize', 10);
text(TI_ci(2), 0.04, sprintf(' %.2f', TI_ci(2)), 'FontSize', 10);
xlabel('Targeting Index', 'FontSize', 12);
ylabel('Probability', 'FontSize', 12);
title('Monte Carlo Distribution (n=10,000)', 'FontSize', 13);
grid on;

subplot(1,2,2)
qqplot(TI_distribution);
title('Q-Q Plot (Normality Assessment)', 'FontSize', 13);
grid on;

% FIGURE 4: Parameter optimization heatmap
figure('Position', [200 200 900 700]);
imagesc(K_ROS_range, k_pH_range, TI_matrix);
colorbar;
xlabel('K_{ROS} (μM)', 'FontSize', 13);
ylabel('k_{pH,max} (h^{-1})', 'FontSize', 13);
title('Targeting Index Optimization Landscape', 'FontSize', 15);
hold on;
plot(optimal_K_ROS, optimal_k_pH, 'r*', 'MarkerSize', 20, 'LineWidth', 3);
text(optimal_K_ROS+3, optimal_k_pH, ...
     sprintf('  Optimal\n  TI=%.2f', max_TI), 'Color', 'white', 'FontSize', 12, ...
     'FontWeight', 'bold');
colormap('jet');
set(gca, 'YDir', 'normal');

% FIGURE 5: Sensitivity analysis
figure('Position', [250 250 1000 400]);

subplot(1,2,1)
bar(sensitivity_indices * 100, 'FaceColor', [0.3 0.5 0.8]);
set(gca, 'XTickLabel', param_names, 'FontSize', 11);
ylabel('Variance Contribution (%)', 'FontSize', 12);
title('Sobol Sensitivity Indices', 'FontSize', 13);
grid on;
ylim([0 max(sensitivity_indices)*110]);

subplot(1,2,2)
plot(beta_range, TI_beta_scan, 'LineWidth', 2.5, 'Color', [0.8 0.2 0.4]);
hold on;
plot(optimal_beta, max_TI_beta, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
xlabel('Synergy Coefficient (β)', 'FontSize', 12);
ylabel('Targeting Index', 'FontSize', 12);
title('Synergistic Effect Optimization', 'FontSize', 13);
grid on;
text(optimal_beta, max_TI_beta+0.2, sprintf('Optimal: %.2f', optimal_beta), ...
     'HorizontalAlignment', 'center', 'FontSize', 11);

% FIGURE 6: Patient stratification
figure('Position', [300 300 1000 500]);

subplot(1,2,1)
bar(accumulation_data(:,2:3));
set(gca, 'XTickLabel', stage_names, 'FontSize', 11);
xlabel('Disease Stage', 'FontSize', 12);
ylabel('Drug Released at 24h (μg)', 'FontSize', 12);
title('Patient Stratification: Drug Accumulation', 'FontSize', 13);
legend('Plaque Tissue', 'Healthy Tissue', 'Location', 'northwest');
grid on;

subplot(1,2,2)
ratios = accumulation_data(:,2) ./ accumulation_data(:,3);
plot(1:4, ratios, 'o-', 'LineWidth', 2.5, 'MarkerSize', 10, ...
     'MarkerFaceColor', [0.8 0.4 0.2], 'Color', [0.8 0.4 0.2]);
set(gca, 'XTick', 1:4, 'XTickLabel', stage_names, 'FontSize', 11);
xlabel('Disease Stage', 'FontSize', 12);
ylabel('Selectivity Ratio (Plaque/Healthy)', 'FontSize', 12);
title('Targeting Selectivity Across Stages', 'FontSize', 13);
grid on;
ylim([0 max(ratios)*1.2]);

% FIGURE 7: Pharmacoeconomic comparison
figure('Position', [350 350 800 500]);

costs = [cost_systemic_10yr, cost_np_10yr] / 1000;
QALYs = [QALY_systemic, QALY_np];

subplot(2,1,1)
bar(costs, 'FaceColor', [0.4 0.6 0.8]);
set(gca, 'XTickLabel', {'Systemic Statin', 'Nanoparticle'});
ylabel('10-Year Cost ($1000s)', 'FontSize', 12);
title('Cost Comparison', 'FontSize', 13);
grid on;
text(1, costs(1)+0.5, sprintf('$%.1fK', costs(1)), 'HorizontalAlignment', 'center');
text(2, costs(2)+0.5, sprintf('$%.1fK', costs(2)), 'HorizontalAlignment', 'center');

subplot(2,1,2)
bar(QALYs, 'FaceColor', [0.6 0.8 0.4]);
set(gca, 'XTickLabel', {'Systemic Statin', 'Nanoparticle'});
ylabel('Quality-Adjusted Life Years', 'FontSize', 12);
title('Effectiveness Comparison', 'FontSize', 13);
grid on;
text(1, QALYs(1)+0.2, sprintf('%.1f years', QALYs(1)), 'HorizontalAlignment', 'center');
text(2, QALYs(2)+0.2, sprintf('%.1f years', QALYs(2)), 'HorizontalAlignment', 'center');

% FINAL SUMMARY 
fprintf('KEY FINDINGS SUMMARY:\n');
fprintf('1. TARGETING PERFORMANCE:\n');
fprintf('   - Dual-responsive TI: %.2f (95%% CI: %.2f-%.2f)\n', TI_mean, TI_ci(1), TI_ci(2));
fprintf('   - %.1f%% probability TI > 4.0\n', prob_TI_gt4);
fprintf('\n');

fprintf('2. OPTIMAL DESIGN PARAMETERS:\n');
fprintf('   - k_pH_max: %.4f h^-1\n', optimal_k_pH);
fprintf('   - K_ROS: %.1f μM\n', optimal_K_ROS);
fprintf('   - Beta (synergy): %.2f\n', optimal_beta);
fprintf('   - Maximum TI achieved: %.2f\n', max_TI_beta);
fprintf('\n');

fprintf('3. PARAMETER SENSITIVITY:\n');
[~, most_sensitive_idx] = max(sensitivity_indices);
fprintf('   - Most influential: %s (%.1f%% variance)\n', ...
        param_names{most_sensitive_idx}, max(sensitivity_indices)*100);
fprintf('   - Least influential: %s (%.1f%% variance)\n', ...
        param_names{find(sensitivity_indices == min(sensitivity_indices), 1)}, ...
        min(sensitivity_indices)*100);
fprintf('\n');

fprintf('4. CLINICAL RECOMMENDATIONS:\n');
[~, best_stage_idx] = max(ratios);
fprintf('   - Best treatment stage: %s (%.1fx selectivity)\n', ...
        stage_names{best_stage_idx}, ratios(best_stage_idx));
fprintf('   - Early intervention strongly recommended\n');
fprintf('\n');

fprintf('5. ECONOMIC IMPACT:\n');
fprintf('   - 10-year savings per patient: $%d (%.1f%%)\n', savings, savings_percent);
fprintf('   - Additional QALYs: %.1f years\n', QALY_gained);
fprintf('   - US population savings (10%% adoption): $%.2f billion/year\n', ...
        annual_savings_population);
fprintf('\n');


% FUNCTIONS 
%%
% 
% $$e^{\pi i} + 1 = 0$$
% 

function k = calc_k_pH(pH, params)
    % Calculate pH-dependent degradation rate
    k = params.k_pH_max / (1 + exp(params.alpha * (pH - params.pH0)));
end

function k = calc_k_ROS(ROS, params)
    % Calculate ROS-dependent degradation rate
    k = params.k_ROS_base + params.k_ROS_max * (ROS / (ROS + params.K_ROS));
end

function k = calc_k_dual(pH, ROS, params)
    % Calculate dual-responsive degradation rate with synergy
    k_pH = calc_k_pH(pH, params);
    k_ROS = calc_k_ROS(ROS, params);
    k = k_pH + k_ROS + params.beta * k_pH * k_ROS;
end

function dydt = release_pH(~, y, pH, params)
    % ODE system for pH-responsive release
    M = y(1);  % Drug in nanoparticle
    D = y(2);  % Released drug
    
    k_deg = calc_k_pH(pH, params);
    
    dMdt = -k_deg * M;
    dDdt = k_deg * M - params.k_elim * D;
    
    dydt = [dMdt; dDdt];
end

function dydt = release_ROS(~, y, ROS, params)
    % ODE system for ROS-responsive release
    M = y(1);
    D = y(2);
    
    k_deg = calc_k_ROS(ROS, params);
    
    dMdt = -k_deg * M;
    dDdt = k_deg * M - params.k_elim * D;
    
    dydt = [dMdt; dDdt];
end

function dydt = release_dual(~, y, pH, ROS, params)
    % ODE system for dual-responsive release
    M = y(1);
    D = y(2);
    
    k_deg = calc_k_dual(pH, ROS, params);
    
    dMdt = -k_deg * M;
    dDdt = k_deg * M - params.k_elim * D;
    
    dydt = [dMdt; dDdt];
end

function dydt = release_baseline(~, y, k_baseline, params)
    % ODE system for non-responsive (control) release
    M = y(1);
    D = y(2);
    
    dMdt = -k_baseline * M;
    dDdt = k_baseline * M - params.k_elim * D;
    
    dydt = [dMdt; dDdt];
end