% ==========================================
% PARAMETRIC STUDY GRAPHS FOR OPTIMIZEDSTUDY2
% ==========================================
% Generates parametric study graphs for the low-pressure-drop design

clear; clc; close all;

fprintf('\n========================================\n');
fprintf('OPTIMIZEDSTUDY2 PARAMETRIC STUDY SUITE\n');
fprintf('Low Pressure Drop Design (ΔP_shell < 1 kPa)\n');
fprintf('========================================\n\n');

% Create output directory
output_dir = '../OptimizedStudy2/Figures/ParametricStudyFigures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Base parameters for OptimizedStudy2 (Design #1)
params.Q = 10000;
params.m_hot = 1.246;
params.m_cold = 1.5;
params.T_hot1 = 293.15;      % 20.00°C EXACTLY
params.T_cold1 = 281.15;     % 8.00°C EXACTLY
params.F = 0.96;
params.L_tube = 2.5;         % OPTIMIZED for Study2
params.D_internal = 0.022;
params.D_external = 0.025;
params.k_tube = 50;
params.l_f = 0.003;          % OPTIMIZED - shorter fins
params.t_f = 0.001;
params.p_f = 0.012;          % OPTIMIZED - wider spacing
params.k_fin = 205;
params.N_p = 4;
params.psi_n = 0.17;
params.C_1 = 0.866;
params.L_bb = 0.0127;
params.baffle_spacing_ratio = 1.5;  % OPTIMIZED - very large
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
params.L_tp = 1.5 * params.D_external;  % OPTIMIZED - wider pitch

% Define studies - ranges centered around OptimizedStudy2 values
studies = {
    {'baffle_spacing_ratio', 1.0, 2.0, 6, 'Baffle Spacing Ratio', '-'};
    {'p_f', 0.008, 0.014, 7, 'Fin Pitch', 'm'};
    {'N_p', 2, 4, 2, 'Number of Tube Passes', '-'};
    {'L_tube', 2.0, 3.5, 6, 'Tube Length', 'm'};
    {'l_f', 0.002, 0.005, 7, 'Fin Height', 'm'};
    {'L_tp', 1.3*params.D_external, 1.6*params.D_external, 4, 'Tube Pitch', 'm'};
};

total_figures = 0;

for study_idx = 1:length(studies)
    study = studies{study_idx};
    var_name = study{1};
    var_min = study{2};
    var_max = study{3};
    num_points = study{4};
    var_label = study{5};
    var_units = study{6};

    fprintf('\n========================================\n');
    fprintf('Study %d/%d: %s\n', study_idx, length(studies), var_label);
    fprintf('========================================\n');

    % Generate test values
    if strcmp(var_name, 'N_p')
        var_values = [2, 4];
    else
        var_values = linspace(var_min, var_max, num_points);
    end

    % Preallocate results
    U_values = zeros(size(var_values));
    Area_values = zeros(size(var_values));
    N_tubes_values = zeros(size(var_values));
    DP_tube_values = zeros(size(var_values));
    DP_shell_values = zeros(size(var_values));
    D_shell_values = zeros(size(var_values));
    converged_flags = zeros(size(var_values));

    % Run calculations
    for i = 1:length(var_values)
        test_params = params;
        test_params.(var_name) = var_values(i);

        % Special handling for L_tp
        if strcmp(var_name, 'L_tp')
            % L_tp is already set correctly
        end

        results = calculateHeatExchanger(test_params);

        U_values(i) = results.U_o_calc;
        Area_values(i) = results.A_o;
        N_tubes_values(i) = ceil(results.N_tt);
        DP_tube_values(i) = results.DP_t_total / 1000;  % Convert to kPa
        DP_shell_values(i) = results.DP_s / 1000;       % Convert to kPa
        D_shell_values(i) = results.D_s;
        converged_flags(i) = results.converged;

        if results.converged
            fprintf('  Point %d/%d: %s = %.4f - Converged\n', i, num_points, var_name, var_values(i));
        else
            fprintf('  Point %d/%d: %s = %.4f - DID NOT CONVERGE\n', i, num_points, var_name, var_values(i));
        end
    end

    % Create figure
    fig = figure('Position', [100, 100, 1200, 800]);

    % Display variable for x-axis
    if strcmp(var_units, 'm')
        x_values_display = var_values * 1000;  % Convert to mm
        x_label = sprintf('%s (mm)', var_label);
    else
        x_values_display = var_values;
        x_label = sprintf('%s (%s)', var_label, var_units);
    end

    % Plot 1: U coefficient
    subplot(2,3,1);
    plot(x_values_display, U_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel(x_label);
    ylabel('U (W/m²K)');
    title('Overall Heat Transfer Coefficient');

    % Plot 2: Area
    subplot(2,3,2);
    plot(x_values_display, Area_values, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel(x_label);
    ylabel('Area (m²)');
    title('Total Heat Transfer Area');

    % Plot 3: Number of tubes
    subplot(2,3,3);
    plot(x_values_display, N_tubes_values, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel(x_label);
    ylabel('Number of Tubes');
    title('Total Number of Tubes');

    % Plot 4: Tube-side pressure drop
    subplot(2,3,4);
    plot(x_values_display, DP_tube_values, 'm-o', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel(x_label);
    ylabel('ΔP_{tube} (kPa)');
    title('Tube-Side Pressure Drop');

    % Plot 5: Shell-side pressure drop (CRITICAL!)
    subplot(2,3,5);
    plot(x_values_display, DP_shell_values, 'k-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    yline(1.0, 'r--', 'LineWidth', 2, 'Label', '1 kPa Limit');
    grid on;
    xlabel(x_label);
    ylabel('ΔP_{shell} (kPa)');
    title('Shell-Side Pressure Drop (CRITICAL)');
    legend('Calculated', '1 kPa Constraint', 'Location', 'best');

    % Plot 6: Shell diameter
    subplot(2,3,6);
    plot(x_values_display, D_shell_values*1000, 'c-o', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel(x_label);
    ylabel('D_{shell} (mm)');
    title('Shell Diameter');

    % Add overall title
    sgtitle(sprintf('OptimizedStudy2 Parametric Study: %s', var_label), 'FontSize', 14, 'FontWeight', 'bold');

    % Save figure
    var_name_clean = strrep(var_name, '_', '');
    filename_base = sprintf('%d_%s', study_idx, var_name_clean);

    saveas(fig, fullfile(output_dir, [filename_base '.png']));
    saveas(fig, fullfile(output_dir, [filename_base '.fig']));

    % Save data
    save(fullfile(output_dir, [filename_base '.mat']), ...
         'var_values', 'U_values', 'Area_values', 'N_tubes_values', ...
         'DP_tube_values', 'DP_shell_values', 'D_shell_values', 'converged_flags');

    fprintf('  Saved figures: %s/%s.png and .fig\n', output_dir, filename_base);
    fprintf('  Saved data: %s/%s.mat\n', output_dir, filename_base);

    total_figures = total_figures + 2;  % PNG + FIG

    close(fig);
end

fprintf('\n========================================\n');
fprintf('ALL STUDIES COMPLETE!\n');
fprintf('========================================\n');
fprintf('Results saved in: %s/\n', output_dir);
fprintf('Total figures generated: %d\n', total_figures);
fprintf('\nFigure files:\n');
fprintf('  - PNG format (for reports/presentations)\n');
fprintf('  - FIG format (editable MATLAB figures)\n');
fprintf('  - MAT format (raw data)\n');
fprintf('========================================\n\n');
