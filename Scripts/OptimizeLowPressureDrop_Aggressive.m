% ==========================================
% AGGRESSIVE OPTIMIZATION FOR ΔP_shell < 1 kPa
% ==========================================
% Previous attempt: minimum 1.42 kPa
% New strategy: MUCH longer tubes, larger spacing, fewer baffles

clear; clc;

% Fixed parameters (CANNOT CHANGE)
params.Q = 10000;
params.m_hot = 1.246;
params.m_cold = 1.5;
params.T_hot1 = 293.15;
params.T_cold1 = 281.15;
params.F = 0.96;
params.D_internal = 0.022;
params.D_external = 0.025;
params.k_tube = 50;
params.k_fin = 205;
params.t_f = 0.001;
params.psi_n = 0.17;
params.C_1 = 0.866;
params.L_bb = 0.0127;
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

fprintf('AGGRESSIVE OPTIMIZATION FOR ΔP_shell < 1 kPa\n');
fprintf('=============================================\n\n');

% More aggressive test matrix
tube_lengths = [2.5, 3.0, 3.5, 4.0];
fin_pitches = [0.008, 0.010, 0.012];  % Larger spacing
tube_pitch_ratios = [1.4, 1.5, 1.6];  % More space
fin_heights = [0.003, 0.004];  % Shorter fins
tube_passes = [2, 4];  % Test fewer passes
baffle_ratios = [1.0, 1.2, 1.5];  % Very large spacing (may mean 0-1 baffles)

best_design = struct();
best_DP = inf;
designs = [];
counter = 0;
total = length(tube_lengths) * length(fin_pitches) * length(tube_pitch_ratios) * ...
        length(fin_heights) * length(tube_passes) * length(baffle_ratios);

fprintf('Testing %d combinations...\n\n', total);

for tl = tube_lengths
    for fp = fin_pitches
        for tpr = tube_pitch_ratios
            for fh = fin_heights
                for np = tube_passes
                    for br = baffle_ratios
                        counter = counter + 1;

                        params.L_tube = tl;
                        params.p_f = fp;
                        params.L_tp = tpr * params.D_external;
                        params.l_f = fh;
                        params.N_p = np;
                        params.baffle_spacing_ratio = br;

                        results = calculateHeatExchanger(params);

                        if results.converged
                            DP_shell_kPa = results.DP_s / 1000;
                            DP_tube_kPa = results.DP_t_total / 1000;

                            design.tube_length = tl;
                            design.fin_pitch = fp;
                            design.tube_pitch_ratio = tpr;
                            design.fin_height = fh;
                            design.N_passes = np;
                            design.baffle_ratio = br;
                            design.U = results.U_o_calc;
                            design.Area = results.A_o;
                            design.N_tubes = ceil(results.N_tt);
                            design.DP_shell = DP_shell_kPa;
                            design.DP_tube = DP_tube_kPa;
                            design.D_s = results.D_s;
                            design.L_B = results.L_B;
                            design.n_baffles = floor(tl / results.L_B) - 1;

                            designs = [designs; design];

                            if DP_shell_kPa < 1.0 && DP_shell_kPa < best_DP
                                best_DP = DP_shell_kPa;
                                best_design = design;
                                fprintf('✓ New best: L=%.1fm, pf=%.0fmm, np=%d, br=%.1f → ΔP_s=%.3f kPa (n_baff=%d)\n', ...
                                        tl, fp*1000, np, br, DP_shell_kPa, design.n_baffles);
                            end
                        end

                        if mod(counter, 50) == 0
                            fprintf('Progress: %d/%d (%.1f%%) | Min ΔP so far: %.3f kPa\n', ...
                                    counter, total, 100*counter/total, min([designs.DP_shell]));
                        end
                    end
                end
            end
        end
    end
end

fprintf('\n=============================================\n');
fprintf('OPTIMIZATION COMPLETE\n');
fprintf('=============================================\n\n');

valid_designs = designs([designs.DP_shell] < 1.0);

if ~isempty(valid_designs)
    fprintf('✓ Found %d designs with ΔP_shell < 1 kPa!\n\n', length(valid_designs));

    [~, idx] = sort([valid_designs.U], 'descend');
    valid_designs = valid_designs(idx);

    fprintf('TOP 10 DESIGNS (sorted by U coefficient):\n');
    fprintf('------------------------------------------------------------------------------------\n');
    fprintf('  # | L_tube | FinPitch | TubePitch | FinHt | Passes | Baffles |    U   | ΔP_s  \n');
    fprintf('    |  (m)   |   (mm)   |   Ratio   | (mm)  |        |         | (W/m²K)| (kPa) \n');
    fprintf('------------------------------------------------------------------------------------\n');

    for i = 1:min(10, length(valid_designs))
        d = valid_designs(i);
        fprintf('%3d |  %.1f   |   %.0f     |   %.2f    | %.1f  |   %d    |    %d    | %6.1f | %.3f\n', ...
                i, d.tube_length, d.fin_pitch*1000, d.tube_pitch_ratio, ...
                d.fin_height*1000, d.N_passes, d.n_baffles, d.U, d.DP_shell);
    end

    fprintf('------------------------------------------------------------------------------------\n\n');

    fprintf('RECOMMENDED DESIGN:\n');
    fprintf('========================================\n');
    d = valid_designs(1);
    fprintf('Tube length: %.1f m\n', d.tube_length);
    fprintf('Fin pitch: %.0f mm (%.0f fins/meter)\n', d.fin_pitch*1000, 1/d.fin_pitch);
    fprintf('Fin height: %.1f mm\n', d.fin_height*1000);
    fprintf('Tube pitch ratio: %.2f (%.1f mm pitch)\n', d.tube_pitch_ratio, d.tube_pitch_ratio*25);
    fprintf('Number of passes: %d\n', d.N_passes);
    fprintf('Baffle spacing ratio: %.2f\n', d.baffle_ratio);
    fprintf('Number of baffles: %d\n', d.n_baffles);
    fprintf('\nPerformance:\n');
    fprintf('  U = %.1f W/m²K\n', d.U);
    fprintf('  Area = %.1f m²\n', d.Area);
    fprintf('  Tubes = %d\n', d.N_tubes);
    fprintf('  Shell diameter = %.0f mm\n', d.D_s*1000);
    fprintf('  ΔP_shell = %.3f kPa ✓✓✓\n', d.DP_shell);
    fprintf('  ΔP_tube = %.2f kPa\n', d.DP_tube);
    fprintf('========================================\n\n');

    best_params = params;
    best_params.L_tube = d.tube_length;
    best_params.p_f = d.fin_pitch;
    best_params.L_tp = d.tube_pitch_ratio * params.D_external;
    best_params.l_f = d.fin_height;
    best_params.N_p = d.N_passes;
    best_params.baffle_spacing_ratio = d.baffle_ratio;

    save('OptimizedStudy2_Design.mat', 'best_params', 'valid_designs');
    fprintf('Saved to: OptimizedStudy2_Design.mat\n');

else
    fprintf('⚠ NO DESIGNS FOUND with ΔP_shell < 1 kPa\n');
    fprintf('Minimum achieved: %.3f kPa\n\n', min([designs.DP_shell]));

    % Show closest designs
    [~, idx] = sort([designs.DP_shell]);
    closest = designs(idx);

    fprintf('CLOSEST 5 DESIGNS:\n');
    fprintf('------------------------------------------------------------------------------------\n');
    for i = 1:min(5, length(closest))
        d = closest(i);
        fprintf('%d. L=%.1fm, pf=%.0fmm, tpr=%.2f, fh=%.1fmm, np=%d, baff=%d → ΔP_s=%.3f kPa\n', ...
                i, d.tube_length, d.fin_pitch*1000, d.tube_pitch_ratio, ...
                d.fin_height*1000, d.N_passes, d.n_baffles, d.DP_shell);
    end
    fprintf('------------------------------------------------------------------------------------\n');
end
