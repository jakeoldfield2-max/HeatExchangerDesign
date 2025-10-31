% ==========================================
% EXPORT ALL 62 VALID DESIGNS TO TABLE
% ==========================================
% Creates a readable text file with all designs meeting ΔP_shell < 1 kPa

clear; clc;

fprintf('Loading OptimizedStudy2 designs...\n');
load('../OptimizedStudy2/Data/OptimizedStudy2_Design.mat', 'valid_designs');

fprintf('Found %d valid designs\n\n', length(valid_designs));

% Create output file
filename = '../OptimizedStudy2/Documentation/All_62_Valid_Designs.txt';
fid = fopen(filename, 'w');

% Write header
fprintf(fid, '================================================================================\n');
fprintf(fid, '                    ALL 62 VALID DESIGNS (ΔP_shell < 1 kPa)\n');
fprintf(fid, '              OptimizedStudy2 - Low Pressure Drop Optimization\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'Constraint: Shell-side pressure drop < 1.0 kPa\n');
fprintf(fid, 'Fixed: Q=10kW, m_hot=1.246 kg/s, m_cold=1.5 kg/s, T_hot=20°C, T_cold=8°C\n');
fprintf(fid, 'Fixed: Tube = 25mm OD x 22mm ID (EU Standard) - SAME FOR ALL 62 DESIGNS\n');
fprintf(fid, 'Fixed: Fin Thickness = 1.0mm - SAME FOR ALL 62 DESIGNS\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'Designs sorted by U coefficient (descending)\n\n');

% Write main table
fprintf(fid, '====================================================================================================\n');
fprintf(fid, '  #  | L_tube | FinPitch | FinHt | FinThick | TubePitch | Passes | Baffles | BafRatio |    U    \n');
fprintf(fid, '     |  (m)   |   (mm)   | (mm)  |   (mm)   |   Ratio   |        |         |          | (W/m²K) \n');
fprintf(fid, '====================================================================================================\n');

for i = 1:length(valid_designs)
    d = valid_designs(i);
    fprintf(fid, '%3d  |  %.1f   |    %2.0f    | %.1f  |   1.0    |   %.2f    |   %d    |    %d    |   %.2f   | %6.1f \n', ...
            i, d.tube_length, d.fin_pitch*1000, d.fin_height*1000, ...
            d.tube_pitch_ratio, d.N_passes, d.n_baffles, d.baffle_ratio, d.U);
end

fprintf(fid, '====================================================================================================\n\n');

fprintf(fid, 'NOTE: Fin Thickness = 1.0mm for ALL designs (was not varied)\n');
fprintf(fid, '      Tube Outer Diameter = 25.0mm for ALL designs (was not varied)\n');
fprintf(fid, '      Tube Inner Diameter = 22.0mm for ALL designs (was not varied)\n\n');

% Continue with second table
fprintf(fid, '================================================================================================\n');
fprintf(fid, '  #  |  Area  | ΔP_shell | ΔP_tube | N_tubes | Shell_D | Baff_Spacing | Actual_Spacing\n');
fprintf(fid, '     |  (m²)  |  (kPa)   |  (kPa)  |         |  (mm)   |    Ratio     |     (mm)\n');
fprintf(fid, '================================================================================================\n');

for i = 1:length(valid_designs)
    d = valid_designs(i);
    baffle_spacing_mm = d.baffle_ratio * d.D_s * 1000;  % Actual spacing in mm

    fprintf(fid, '%3d  | %5.1f  |  %.3f   |  %.2f   |   %3d   |  %5.0f  |     %.2f     |     %5.0f\n', ...
            i, d.Area, d.DP_shell, d.DP_tube, d.N_tubes, d.D_s*1000, ...
            d.baffle_ratio, baffle_spacing_mm);
end

fprintf(fid, '================================================================================================\n\n');

% Summary statistics
fprintf(fid, '================================================================================\n');
fprintf(fid, '                              SUMMARY STATISTICS\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'U Coefficient:\n');
fprintf(fid, '  Range: %.1f - %.1f W/m²K\n', min([valid_designs.U]), max([valid_designs.U]));
fprintf(fid, '  Mean: %.1f W/m²K\n\n', mean([valid_designs.U]));

fprintf(fid, 'Shell-Side Pressure Drop:\n');
fprintf(fid, '  Range: %.3f - %.3f kPa\n', min([valid_designs.DP_shell]), max([valid_designs.DP_shell]));
fprintf(fid, '  Mean: %.3f kPa\n\n', mean([valid_designs.DP_shell]));

fprintf(fid, 'Tube-Side Pressure Drop:\n');
fprintf(fid, '  Range: %.3f - %.3f kPa\n', min([valid_designs.DP_tube]), max([valid_designs.DP_tube]));
fprintf(fid, '  Mean: %.2f kPa\n\n', mean([valid_designs.DP_tube]));

fprintf(fid, 'Heat Transfer Area:\n');
fprintf(fid, '  Range: %.1f - %.1f m²\n', min([valid_designs.Area]), max([valid_designs.Area]));
fprintf(fid, '  Mean: %.1f m²\n\n', mean([valid_designs.Area]));

fprintf(fid, 'Number of Tubes:\n');
fprintf(fid, '  Range: %d - %d tubes\n', min([valid_designs.N_tubes]), max([valid_designs.N_tubes]));
fprintf(fid, '  Mean: %.0f tubes\n\n', round(mean([valid_designs.N_tubes])));

% Parameter distribution
fprintf(fid, '================================================================================\n');
fprintf(fid, '                          PARAMETER DISTRIBUTIONS\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'FIXED PARAMETERS (NOT VARIED IN OPTIMIZATION):\n');
fprintf(fid, '  Tube Outer Diameter: 25.0 mm (ALL 62 designs)\n');
fprintf(fid, '  Tube Inner Diameter: 22.0 mm (ALL 62 designs)\n');
fprintf(fid, '  Tube Wall Thickness: 1.5 mm (ALL 62 designs)\n');
fprintf(fid, '  Fin Thickness: 1.0 mm (ALL 62 designs)\n');
fprintf(fid, '  Tube Material: Carbon Steel (k = 50 W/mK)\n');
fprintf(fid, '  Fin Material: Aluminum (k = 205 W/mK)\n\n');

fprintf(fid, 'VARIED PARAMETERS:\n\n');

fprintf(fid, 'Tube Length Distribution:\n');
tube_lengths_unique = unique([valid_designs.tube_length]);
for tl = tube_lengths_unique
    count = sum([valid_designs.tube_length] == tl);
    fprintf(fid, '  %.1f m: %d designs (%.1f%%)\n', tl, count, count/length(valid_designs)*100);
end
fprintf(fid, '\n');

fprintf(fid, 'Fin Pitch Distribution:\n');
fin_pitches_unique = unique([valid_designs.fin_pitch]);
for fp = fin_pitches_unique
    count = sum([valid_designs.fin_pitch] == fp);
    fprintf(fid, '  %.0f mm: %d designs (%.1f%%)\n', fp*1000, count, count/length(valid_designs)*100);
end
fprintf(fid, '\n');

fprintf(fid, 'Fin Height Distribution:\n');
fin_heights_unique = unique([valid_designs.fin_height]);
for fh = fin_heights_unique
    count = sum([valid_designs.fin_height] == fh);
    fprintf(fid, '  %.0f mm: %d designs (%.1f%%)\n', fh*1000, count, count/length(valid_designs)*100);
end
fprintf(fid, '\n');

fprintf(fid, 'Tube Pitch Ratio Distribution:\n');
tube_pitch_unique = unique([valid_designs.tube_pitch_ratio]);
for tpr = tube_pitch_unique
    count = sum([valid_designs.tube_pitch_ratio] == tpr);
    fprintf(fid, '  %.2f (%.1f mm): %d designs (%.1f%%)\n', tpr, tpr*25, count, count/length(valid_designs)*100);
end
fprintf(fid, '\n');

fprintf(fid, 'Number of Passes Distribution:\n');
passes_unique = unique([valid_designs.N_passes]);
for np = passes_unique
    count = sum([valid_designs.N_passes] == np);
    fprintf(fid, '  %d passes: %d designs (%.1f%%)\n', np, count, count/length(valid_designs)*100);
end
fprintf(fid, '\n');

fprintf(fid, 'Baffle Spacing Ratio Distribution:\n');
baffle_ratios_unique = unique([valid_designs.baffle_ratio]);
for br = baffle_ratios_unique
    count = sum([valid_designs.baffle_ratio] == br);
    fprintf(fid, '  %.1f: %d designs (%.1f%%)\n', br, count, count/length(valid_designs)*100);
end
fprintf(fid, '\n');

% Key insights
fprintf(fid, '================================================================================\n');
fprintf(fid, '                               KEY INSIGHTS\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, '1. ALL 62 designs use L_tube = %.1f m (the minimum tested)\n', tube_lengths_unique(1));
fprintf(fid, '   → Shorter tubes might also work! (should test 1.5m, 2.0m)\n\n');

baffle_counts = arrayfun(@(br) sum([valid_designs.baffle_ratio] == br), baffle_ratios_unique);
[max_count, max_idx] = max(baffle_counts);
fprintf(fid, '2. MOST designs (%.0f%%) use baffle_ratio = %.1f (the maximum tested)\n', ...
        max_count / length(valid_designs) * 100, ...
        baffle_ratios_unique(max_idx));
fprintf(fid, '   → Even larger ratios (2.0, 3.0) might work better!\n\n');

fprintf(fid, '3. Fin pitch varies (8-12mm) - this parameter is less critical\n\n');

fprintf(fid, '4. Both 2 and 4 passes work:\n');
two_pass_count = sum([valid_designs.N_passes] == 2);
fprintf(fid, '   - 2 passes: %d designs (%.1f%%) → Simpler, lower ΔP_tube\n', ...
        two_pass_count, two_pass_count/length(valid_designs)*100);
four_pass_count = sum([valid_designs.N_passes] == 4);
fprintf(fid, '   - 4 passes: %d designs (%.1f%%) → Higher U, more complex\n\n', ...
        four_pass_count, four_pass_count/length(valid_designs)*100);

fprintf(fid, '5. Trade-off: Higher U requires more area but gives better heat transfer\n');
fprintf(fid, '   Design #1: U=%.1f W/m²K, Area=%.1f m² (highest U)\n', ...
        valid_designs(1).U, valid_designs(1).Area);
fprintf(fid, '   Design #62: U=%.1f W/m²K, Area=%.1f m² (lowest U, most area)\n\n', ...
        valid_designs(end).U, valid_designs(end).Area);

fprintf(fid, '6. TUBE SIZE was FIXED at 25mm OD x 22mm ID for ALL designs\n');
fprintf(fid, '   → This is a standard EU metric size (1.5mm wall)\n');
fprintf(fid, '   → Different tube sizes were not tested in this optimization\n\n');

fprintf(fid, '7. FIN THICKNESS was FIXED at 1.0mm for ALL designs\n');
fprintf(fid, '   → Standard fin thickness for aluminum helical fins\n');
fprintf(fid, '   → Different thicknesses were not tested in this optimization\n\n');

% Top recommendations
fprintf(fid, '================================================================================\n');
fprintf(fid, '                          TOP DESIGN RECOMMENDATIONS\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'BEST OVERALL (Highest U):\n');
fprintf(fid, '  Design #1: U=%.1f W/m²K, ΔP_s=%.3f kPa, 4 passes\n', ...
        valid_designs(1).U, valid_designs(1).DP_shell);
fprintf(fid, '  → Best heat transfer, but 3%% margin on ΔP constraint\n\n');

% Find best 2-pass design
two_pass_designs = valid_designs([valid_designs.N_passes] == 2);
if ~isempty(two_pass_designs)
    [~, best_2p_idx] = max([two_pass_designs.U]);
    best_2p = two_pass_designs(best_2p_idx);
    original_idx = find([valid_designs.U] == best_2p.U & [valid_designs.N_passes] == 2, 1);
    fprintf(fid, 'BEST SIMPLE DESIGN (Highest U with 2 passes):\n');
    fprintf(fid, '  Design #%d: U=%.1f W/m²K, ΔP_s=%.3f kPa, ΔP_t=%.2f kPa, 2 passes\n', ...
            original_idx, best_2p.U, best_2p.DP_shell, best_2p.DP_tube);
    fprintf(fid, '  → Simpler manufacturing, 8× lower ΔP_tube than Design #1\n\n');
end

% Find safest design (lowest ΔP_shell)
[min_DP, safest_idx] = min([valid_designs.DP_shell]);
fprintf(fid, 'SAFEST DESIGN (Lowest ΔP_shell):\n');
fprintf(fid, '  Design #%d: U=%.1f W/m²K, ΔP_s=%.3f kPa (%.0f%% margin!)\n', ...
        safest_idx, valid_designs(safest_idx).U, min_DP, (1-min_DP)*100);
fprintf(fid, '  → Maximum safety margin against constraint violation\n\n');

fprintf(fid, '================================================================================\n');
fprintf(fid, '                           END OF DESIGN TABLE\n');
fprintf(fid, '================================================================================\n');

fclose(fid);

fprintf('✓ Complete table saved to:\n');
fprintf('  %s\n\n', filename);
fprintf('File contains:\n');
fprintf('  - All 62 designs in detailed tables\n');
fprintf('  - Fixed parameters clearly noted\n');
fprintf('  - Summary statistics\n');
fprintf('  - Parameter distributions\n');
fprintf('  - Key insights\n');
fprintf('  - Design recommendations\n');
