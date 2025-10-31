% ==========================================
% RUN ALL PARAMETRIC STUDIES AND SAVE GRAPHS
% ==========================================
% Generates comprehensive parametric study results with saved figures

clear; clc; close all;

fprintf('\n========================================\n');
fprintf('COMPREHENSIVE PARAMETRIC STUDY SUITE\n');
fprintf('========================================\n\n');

% Create output directory for figures
if ~exist('ParametricStudyFigures', 'dir')
    mkdir('ParametricStudyFigures');
end

% Base parameters (OPTIMIZED with CORRECTED TEMPERATURES)
params.Q = 10000;
params.m_hot = 1.246;
params.m_cold = 1.5;
params.T_hot1 = 293.15;      % 20.00°C EXACTLY
params.T_cold1 = 281.15;     % 8.00°C EXACTLY
params.F = 0.96;
params.L_tube = 1.0;
params.D_internal = 0.022;
params.D_external = 0.025;
params.k_tube = 50;
params.l_f = 0.005;
params.t_f = 0.001;
params.p_f = 0.005;
params.k_fin = 205;
params.N_p = 4;
params.psi_n = 0.17;
params.C_1 = 0.866;
params.L_bb = 0.0127;
params.baffle_spacing_ratio = 0.7;
params.Rf_i = 0.00018;
params.Rf_o = 0.0003526;
params.U_initial_guess = 80;
params.tolerance = 0.1;
params.j_f = 0.04;
params.rho_t = 999.7;
params.mu_t = 1.308e-3;
params.k_t = 0.58;
params.Cp_t = 4200;
params.rho_s = 1.204;
params.mu_s = 1.825e-5;
params.k_s = 0.02514;
params.Cp_s = 1007;
params.L_tp = 1.25 * params.D_external;

% Define all parametric studies
studies = {
    % {variable_name, min, max, num_points, display_name, units}
    {'baffle_spacing_ratio', 0.3, 0.8, 6, 'Baffle Spacing Ratio', '-'};
    {'p_f', 0.003, 0.006, 7, 'Fin Pitch', 'm'};
    {'N_p', 2, 8, 4, 'Number of Tube Passes', '-'};
    {'D_internal', 0.019, 0.022, 4, 'Tube Inner Diameter', 'm'};
    {'L_tube', 1.0, 2.0, 5, 'Tube Length', 'm'};
    {'l_f', 0.003, 0.008, 6, 'Fin Height', 'm'};
    {'t_f', 0.0005, 0.002, 6, 'Fin Thickness', 'm'};
};

% Run each study
for study_idx = 1:length(studies)
    study = studies{study_idx};
    var_name = study{1};
    min_val = study{2};
    max_val = study{3};
    num_pts = study{4};
    display_name = study{5};
    units = study{6};

    fprintf('\n========================================\n');
    fprintf('Study %d/%d: %s\n', study_idx, length(studies), display_name);
    fprintf('========================================\n');

    % Generate test values
    test_values = linspace(min_val, max_val, num_pts);

    % Preallocate results
    U_vals = zeros(num_pts, 1);
    Area_vals = zeros(num_pts, 1);
    N_tt_vals = zeros(num_pts, 1);
    DP_t_vals = zeros(num_pts, 1);
    DP_s_vals = zeros(num_pts, 1);
    D_s_vals = zeros(num_pts, 1);
    converged_vals = zeros(num_pts, 1);

    % Run calculations
    for i = 1:num_pts
        test_params = params;
        test_params.(var_name) = test_values(i);

        % Special cases
        if strcmp(var_name, 'D_external')
            test_params.L_tp = 1.25 * test_params.D_external;
        end

        try
            results = calculateHeatExchanger(test_params);
            U_vals(i) = results.U_o_calc;
            Area_vals(i) = results.A_o;
            N_tt_vals(i) = results.N_tt;
            DP_t_vals(i) = results.DP_t_total / 1000;
            DP_s_vals(i) = results.DP_s / 1000;
            D_s_vals(i) = results.D_s;
            converged_vals(i) = results.converged;
        catch
            U_vals(i) = NaN;
            Area_vals(i) = NaN;
            N_tt_vals(i) = NaN;
            DP_t_vals(i) = NaN;
            DP_s_vals(i) = NaN;
            D_s_vals(i) = NaN;
            converged_vals(i) = 0;
        end
    end

    fprintf('  Completed %d calculations\n', num_pts);

    % Create figure with 6 subplots
    fig = figure('Position', [100, 100, 1400, 900]);

    % Subplot 1: U coefficient
    subplot(2, 3, 1);
    plot(test_values, U_vals, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel(sprintf('%s (%s)', display_name, units));
    ylabel('U (W/m²K)');
    title('Overall Heat Transfer Coefficient');
    grid on;

    % Subplot 2: Total Area
    subplot(2, 3, 2);
    plot(test_values, Area_vals, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel(sprintf('%s (%s)', display_name, units));
    ylabel('Area (m²)');
    title('Total Heat Transfer Area');
    grid on;

    % Subplot 3: Number of Tubes
    subplot(2, 3, 3);
    plot(test_values, ceil(N_tt_vals), '-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel(sprintf('%s (%s)', display_name, units));
    ylabel('Number of Tubes');
    title('Total Number of Tubes');
    grid on;

    % Subplot 4: Tube-Side Pressure Drop
    subplot(2, 3, 4);
    plot(test_values, DP_t_vals, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel(sprintf('%s (%s)', display_name, units));
    ylabel('Pressure Drop (kPa)');
    title('Tube-Side Pressure Drop');
    grid on;

    % Subplot 5: Shell-Side Pressure Drop
    subplot(2, 3, 5);
    plot(test_values, DP_s_vals, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel(sprintf('%s (%s)', display_name, units));
    ylabel('Pressure Drop (kPa)');
    title('Shell-Side Pressure Drop');
    grid on;

    % Subplot 6: Shell Diameter
    subplot(2, 3, 6);
    plot(test_values, D_s_vals, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel(sprintf('%s (%s)', display_name, units));
    ylabel('Shell Diameter (m)');
    title('Shell Diameter');
    grid on;

    % Add overall title
    sgtitle(sprintf('Parametric Study: Effect of %s', display_name), 'FontSize', 14, 'FontWeight', 'bold');

    % Save figure
    filename = sprintf('ParametricStudyFigures/%d_%s', study_idx, strrep(var_name, '_', ''));
    saveas(fig, [filename '.png']);
    saveas(fig, [filename '.fig']);

    fprintf('  Saved figures: %s.png and %s.fig\n', filename, filename);

    % Save data
    save([filename '.mat'], 'test_values', 'U_vals', 'Area_vals', 'N_tt_vals', ...
         'DP_t_vals', 'DP_s_vals', 'D_s_vals', 'converged_vals', 'var_name', 'display_name');

    fprintf('  Saved data: %s.mat\n', filename);
end

fprintf('\n========================================\n');
fprintf('ALL STUDIES COMPLETE!\n');
fprintf('========================================\n');
fprintf('Results saved in: ParametricStudyFigures/\n');
fprintf('Total figures generated: %d\n', length(studies) * 2);
fprintf('\nFigure files:\n');
fprintf('  - PNG format (for reports/presentations)\n');
fprintf('  - FIG format (editable MATLAB figures)\n');
fprintf('  - MAT format (raw data)\n');
fprintf('========================================\n\n');
