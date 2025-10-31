% ==========================================
% === HEAT EXCHANGER DESIGN SPECIFICATION ==
% ==========================================
%
% This script generates a comprehensive design specification document
% for your heat exchanger design. It creates a text file with:
%   - Key performance metrics
%   - All design parameters
%   - Physical dimensions
%   - Operating conditions
%   - Design implications
%
% HOW TO USE:
% 1. Modify the parameters in this script (or it will use defaults)
% 2. Run the script
% 3. A file "HeatExchanger_DesignSpec_[timestamp].txt" will be created
%
% NOTE: You can also copy parameters from ParametricStudy.m

clear; clc;

% ==========================================
% ====== DEFINE DESIGN PARAMETERS ==========
% ==========================================

% --- Specifications
params.Q = 10000;           % Heat duty (W)
params.m_hot = 1.246;       % Air mass flow (kg/s)
params.m_cold = 1.5;        % Water mass flow (kg/s)
params.T_hot1 = 293.15;     % Air in (K) - EXACTLY 20.00°C
params.T_cold1 = 281.15;    % Water in (K) - EXACTLY 8.00°C
params.F = 0.96;            % LMTD Correction factor

% --- Geometry (OPTIMIZED - EU METRIC STANDARD)
params.L_tube = 1.0;        % Tube length (m) - OPTIMIZED
params.D_internal = 0.022;  % Tube inner diameter (m) - 22mm STANDARD
params.D_external = 0.025;  % Tube outer diameter (m) - 25mm STANDARD
params.k_tube = 50;         % Tube wall conductivity (W/mK)

% --- Fins (OPTIMIZED)
params.l_f = 0.005;         % Fin height (m)
params.t_f = 0.001;         % Fin thickness (m)
params.p_f = 0.005;         % Fin pitch (m) - OPTIMIZED (200 fins/meter)
params.k_fin = 205;         % Fin conductivity (W/mK)

% --- Layout (OPTIMIZED)
params.N_p = 4;             % Number of tube passes - OPTIMIZED
params.psi_n = 0.17;        % Layout parameter
params.C_1 = 0.866;         % Layout constant (triangular 30deg)
params.L_bb = 0.0127;       % Baffle bypass clearance (m)
params.baffle_spacing_ratio = 0.8;  % Baffle spacing - OPTIMIZED

% --- Fouling
params.Rf_i = 0.00018;      % Tube-side fouling (m2K/W)
params.Rf_o = 0.0003526;    % Shell-side fouling (m2K/W)

% --- Calculation Parameters
params.U_initial_guess = 80;  % Initial guess for U (W/m2K)
params.tolerance = 0.1;       % Convergence tolerance (W/m2K)
params.j_f = 0.04;            % Shell-side friction factor (from chart)

% --- Tube Side Fluid Properties (Water @ 10C)
params.rho_t = 999.7;       % Density (kg/m3)
params.mu_t = 1.308e-3;     % Viscosity (Pa*s)
params.k_t = 0.58;          % Conductivity (W/mK)
params.Cp_t = 4200;         % Spec. Heat (J/kgK)

% --- Shell Side Fluid Properties (Air @ 20C)
params.rho_s = 1.204;       % Density (kg/m3)
params.mu_s = 1.825e-5;     % Viscosity (Pa*s)
params.k_s = 0.02514;       % Conductivity (W/mK)
params.Cp_s = 1007;         % Spec. Heat (J/kgK)

% Calculate L_tp based on D_external
params.L_tp = 1.25 * params.D_external;

% ==========================================
% ======= RUN HEAT EXCHANGER DESIGN ========
% ==========================================

fprintf('Running heat exchanger calculation...\n');
results = calculateHeatExchanger(params);

if ~results.converged
    warning('Design did not converge! Results may be inaccurate.');
end

% Calculate additional useful parameters
T_hot2 = params.T_hot1 - params.Q / (params.m_hot * params.Cp_s);
T_cold2 = params.T_cold1 + params.Q / (params.m_cold * params.Cp_t);
LMTD = ((params.T_hot1 - T_cold2) - (T_hot2 - params.T_cold1)) / ...
       log((params.T_hot1 - T_cold2)/(T_hot2 - params.T_cold1));

% ==========================================
% ===== GENERATE SPECIFICATION DOCUMENT ====
% ==========================================

% Create filename with timestamp
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
filename = sprintf('HeatExchanger_DesignSpec_%s.txt', timestamp);

fprintf('Generating design specification document...\n');

% Open file for writing
fid = fopen(filename, 'w');

% Write header
fprintf(fid, '================================================================================\n');
fprintf(fid, '              HEAT EXCHANGER DESIGN SPECIFICATION DOCUMENT\n');
fprintf(fid, '================================================================================\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Design Type: Shell-and-Tube Heat Exchanger with Helical Fins\n');
fprintf(fid, 'Application: Air-to-Water Heat Exchange\n');
fprintf(fid, 'Convergence: %s (%d iterations)\n', ...
        iif(results.converged, 'SUCCESS', 'FAILED'), results.iterations);
fprintf(fid, '================================================================================\n\n');

% ========================================
% KEY PERFORMANCE METRICS
% ========================================
fprintf(fid, '################################################################################\n');
fprintf(fid, '#                         KEY PERFORMANCE METRICS                              #\n');
fprintf(fid, '################################################################################\n\n');

fprintf(fid, '--- THERMAL PERFORMANCE ---\n');
fprintf(fid, '  Overall Heat Transfer Coefficient (U):  %.2f W/m²K\n', results.U_o_calc);
fprintf(fid, '  Total Heat Transfer Area:                %.2f m²\n', results.A_o);
fprintf(fid, '  Heat Duty:                               %.2f kW\n', params.Q/1000);
fprintf(fid, '  Log Mean Temperature Difference (LMTD):  %.2f K\n', LMTD);
fprintf(fid, '  Effectiveness (ε):                       %.3f\n', ...
        params.Q/(min(params.m_hot*params.Cp_s, params.m_cold*params.Cp_t)*(params.T_hot1-params.T_cold1)));
fprintf(fid, '\n');

fprintf(fid, '--- PHYSICAL SIZE ---\n');
fprintf(fid, '  Number of Tubes:                         %d\n', ceil(results.N_tt));
fprintf(fid, '  Tube Length:                             %.3f m\n', params.L_tube);
fprintf(fid, '  Shell Diameter:                          %.3f m (%.1f mm)\n', results.D_s, results.D_s*1000);
fprintf(fid, '  Baffle Spacing:                          %.3f m (%.1f mm)\n', results.L_B, results.L_B*1000);
fprintf(fid, '  Overall Length:                          ~%.3f m\n', params.L_tube + 0.3);
fprintf(fid, '\n');

fprintf(fid, '--- HYDRAULIC PERFORMANCE ---\n');
fprintf(fid, '  Tube-Side Pressure Drop:                 %.2f kPa\n', results.DP_t_total/1000);
fprintf(fid, '  Shell-Side Pressure Drop:                %.2f kPa\n', results.DP_s/1000);
fprintf(fid, '  Tube-Side Reynolds Number:               %.0f (%s)\n', results.Re_t, ...
        iif(results.Re_t > 10000, 'Turbulent', iif(results.Re_t > 2300, 'Transitional', 'Laminar')));
fprintf(fid, '  Shell-Side Reynolds Number:              %.0f (%s)\n', results.Re_s, ...
        iif(results.Re_s > 10000, 'Turbulent', iif(results.Re_s > 2300, 'Transitional', 'Laminar')));
fprintf(fid, '\n');

fprintf(fid, '--- TEMPERATURES ---\n');
fprintf(fid, '  Hot Fluid Inlet:                         %.2f K (%.2f °C)\n', ...
        params.T_hot1, params.T_hot1-273.15);
fprintf(fid, '  Hot Fluid Outlet:                        %.2f K (%.2f °C)\n', ...
        T_hot2, T_hot2-273.15);
fprintf(fid, '  Hot Fluid Temperature Drop:              %.2f K\n', params.T_hot1 - T_hot2);
fprintf(fid, '  Cold Fluid Inlet:                        %.2f K (%.2f °C)\n', ...
        params.T_cold1, params.T_cold1-273.15);
fprintf(fid, '  Cold Fluid Outlet:                       %.2f K (%.2f °C)\n', ...
        T_cold2, T_cold2-273.15);
fprintf(fid, '  Cold Fluid Temperature Rise:             %.2f K\n', T_cold2 - params.T_cold1);
fprintf(fid, '\n\n');

% ========================================
% DETAILED DESIGN PARAMETERS
% ========================================
fprintf(fid, '################################################################################\n');
fprintf(fid, '#                       DETAILED DESIGN PARAMETERS                             #\n');
fprintf(fid, '################################################################################\n\n');

% --- Operating Conditions
fprintf(fid, '========================================\n');
fprintf(fid, 'OPERATING CONDITIONS\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Heat Duty (Q):                    %.2f kW\n', params.Q/1000);
fprintf(fid, '  → Implication: Total thermal power to be transferred\n\n');
fprintf(fid, 'Hot Side Mass Flow (m_hot):       %.3f kg/s\n', params.m_hot);
fprintf(fid, '  → Implication: Higher flow → better heat transfer, higher pressure drop\n\n');
fprintf(fid, 'Cold Side Mass Flow (m_cold):     %.3f kg/s\n', params.m_cold);
fprintf(fid, '  → Implication: Higher flow → better heat transfer, higher pressure drop\n\n');
fprintf(fid, 'LMTD Correction Factor (F):       %.3f\n', params.F);
fprintf(fid, '  → Implication: Accounts for multi-pass configuration deviation from counterflow\n');
fprintf(fid, '  → F = %.3f indicates %s configuration effectiveness\n', params.F, ...
        iif(params.F > 0.9, 'excellent', iif(params.F > 0.8, 'good', 'acceptable')));
fprintf(fid, '\n\n');

% --- Tube Geometry
fprintf(fid, '========================================\n');
fprintf(fid, 'TUBE GEOMETRY\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Tube Length (L_tube):             %.3f m\n', params.L_tube);
fprintf(fid, '  → Implication: Longer tubes → more area but higher pressure drop\n\n');
fprintf(fid, 'Tube Inner Diameter (D_internal): %.3f m (%.1f mm)\n', params.D_internal, params.D_internal*1000);
fprintf(fid, '  → Implication: Larger bore → lower velocity, lower pressure drop, lower h_i\n\n');
fprintf(fid, 'Tube Outer Diameter (D_external): %.3f m (%.1f mm)\n', params.D_external, params.D_external*1000);
fprintf(fid, '  → Implication: Determines fin base diameter and shell-side flow area\n\n');
fprintf(fid, 'Wall Thickness:                   %.4f m (%.2f mm)\n', ...
        (params.D_external-params.D_internal)/2, (params.D_external-params.D_internal)/2*1000);
fprintf(fid, '  → Implication: Provides structural integrity and conduction path\n\n');
fprintf(fid, 'Tube Material Conductivity:       %.1f W/mK\n', params.k_tube);
fprintf(fid, '  → Material: %s\n', iif(params.k_tube > 300, 'Copper', iif(params.k_tube > 150, 'Aluminum', 'Steel')));
fprintf(fid, '  → Implication: Higher conductivity → lower wall thermal resistance\n\n\n');

% --- Fin Parameters
fprintf(fid, '========================================\n');
fprintf(fid, 'FIN PARAMETERS (Helical Fins)\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Fin Height (l_f):                 %.3f m (%.1f mm)\n', params.l_f, params.l_f*1000);
fprintf(fid, '  → Implication: Taller fins increase area but may reduce fin efficiency\n');
fprintf(fid, '  → Fin efficiency depends on h_s, k_fin, and fin geometry\n\n');
fprintf(fid, 'Fin Thickness (t_f):              %.4f m (%.2f mm)\n', params.t_f, params.t_f*1000);
fprintf(fid, '  → Implication: Thicker fins → better efficiency but fewer fins per meter\n\n');
fprintf(fid, 'Fin Pitch (p_f):                  %.4f m (%.2f mm)\n', params.p_f, params.p_f*1000);
fprintf(fid, '  → Fins per meter: %.1f\n', 1/params.p_f);
fprintf(fid, '  → Total fins per tube: %.0f\n', params.L_tube/params.p_f);
fprintf(fid, '  → Implication: Smaller pitch → more fins, more area, higher shell-side ΔP\n\n');
fprintf(fid, 'Fin Material Conductivity:        %.1f W/mK\n', params.k_fin);
fprintf(fid, '  → Material: %s\n', iif(params.k_fin > 300, 'Copper', iif(params.k_fin > 150, 'Aluminum', 'Steel')));
fprintf(fid, '  → Implication: Higher conductivity → better fin efficiency\n\n\n');

% --- Layout Configuration
fprintf(fid, '========================================\n');
fprintf(fid, 'LAYOUT & CONFIGURATION\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Number of Tube Passes (N_p):      %d\n', params.N_p);
fprintf(fid, '  → Tubes per pass: %.1f\n', results.N_tt/params.N_p);
fprintf(fid, '  → Total flow path length: %.2f m\n', params.L_tube * params.N_p);
fprintf(fid, '  → Implication: More passes → higher velocity → better h_i but higher ΔP_t\n\n');
fprintf(fid, 'Tube Pitch (L_tp):                %.4f m (%.2f mm)\n', params.L_tp, params.L_tp*1000);
fprintf(fid, '  → Pitch/Diameter ratio: %.2f\n', params.L_tp/params.D_external);
fprintf(fid, '  → Implication: Larger pitch → less compact bundle, lower shell-side ΔP\n\n');
fprintf(fid, 'Tube Arrangement:                 %s\n', iif(params.C_1 == 0.866, 'Triangular (30°)', 'Square'));
fprintf(fid, '  → Layout constant (C_1): %.3f\n', params.C_1);
fprintf(fid, '  → Implication: Triangular is more compact, square easier to clean\n\n');
fprintf(fid, 'Baffle Spacing Ratio:             %.2f\n', params.baffle_spacing_ratio);
fprintf(fid, '  → Actual baffle spacing: %.3f m (%.1f mm)\n', results.L_B, results.L_B*1000);
fprintf(fid, '  → Number of baffles: ~%d\n', floor(params.L_tube/results.L_B) - 1);
fprintf(fid, '  → Implication: Smaller spacing → more baffles → better h_s but much higher ΔP_s\n\n');
fprintf(fid, 'Baffle Bypass Clearance (L_bb):   %.4f m (%.2f mm)\n', params.L_bb, params.L_bb*1000);
fprintf(fid, '  → Implication: Clearance between tube bundle and shell for bypass flow\n\n\n');

% --- Fouling Allowances
fprintf(fid, '========================================\n');
fprintf(fid, 'FOULING ALLOWANCES\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Tube-Side Fouling (Rf_i):         %.6f m²K/W\n', params.Rf_i);
fprintf(fid, '  → Fluid: Water\n');
fprintf(fid, '  → Implication: Design margin for deposit buildup over time\n\n');
fprintf(fid, 'Shell-Side Fouling (Rf_o):        %.6f m²K/W\n', params.Rf_o);
fprintf(fid, '  → Fluid: Air\n');
fprintf(fid, '  → Implication: Design margin for particulate/dust accumulation\n\n\n');

% --- Fluid Properties
fprintf(fid, '========================================\n');
fprintf(fid, 'FLUID PROPERTIES\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'TUBE SIDE (Water at ~%.0f°C):\n', params.T_cold1-273.15);
fprintf(fid, '  Density (rho_t):       %.2f kg/m³\n', params.rho_t);
fprintf(fid, '  Viscosity (mu_t):      %.3e Pa·s\n', params.mu_t);
fprintf(fid, '  Conductivity (k_t):    %.3f W/mK\n', params.k_t);
fprintf(fid, '  Specific Heat (Cp_t):  %.0f J/kgK\n', params.Cp_t);
fprintf(fid, '  Prandtl Number (Pr_t): %.2f\n', params.Cp_t*params.mu_t/params.k_t);
fprintf(fid, '\n');
fprintf(fid, 'SHELL SIDE (Air at ~%.0f°C):\n', params.T_hot1-273.15);
fprintf(fid, '  Density (rho_s):       %.3f kg/m³\n', params.rho_s);
fprintf(fid, '  Viscosity (mu_s):      %.3e Pa·s\n', params.mu_s);
fprintf(fid, '  Conductivity (k_s):    %.5f W/mK\n', params.k_s);
fprintf(fid, '  Specific Heat (Cp_s):  %.0f J/kgK\n', params.Cp_s);
fprintf(fid, '  Prandtl Number (Pr_s): %.3f\n', params.Cp_s*params.mu_s/params.k_s);
fprintf(fid, '\n\n');

% --- Design Notes
fprintf(fid, '################################################################################\n');
fprintf(fid, '#                            DESIGN NOTES                                      #\n');
fprintf(fid, '################################################################################\n\n');

fprintf(fid, '1. PRESSURE DROP ASSESSMENT:\n');
if results.DP_t_total/1000 < 20
    fprintf(fid, '   Tube-side ΔP (%.2f kPa): ACCEPTABLE - Low pumping cost\n', results.DP_t_total/1000);
elseif results.DP_t_total/1000 < 50
    fprintf(fid, '   Tube-side ΔP (%.2f kPa): MODERATE - Acceptable for most applications\n', results.DP_t_total/1000);
else
    fprintf(fid, '   Tube-side ΔP (%.2f kPa): HIGH - Consider increasing D_internal or reducing N_p\n', results.DP_t_total/1000);
end

if results.DP_s/1000 < 100
    fprintf(fid, '   Shell-side ΔP (%.2f kPa): ACCEPTABLE - Low blower power\n', results.DP_s/1000);
elseif results.DP_s/1000 < 500
    fprintf(fid, '   Shell-side ΔP (%.2f kPa): MODERATE - Acceptable for most applications\n', results.DP_s/1000);
else
    fprintf(fid, '   Shell-side ΔP (%.2f kPa): HIGH - Consider increasing baffle spacing or p_f\n', results.DP_s/1000);
end
fprintf(fid, '\n');

fprintf(fid, '2. FLOW REGIME ASSESSMENT:\n');
fprintf(fid, '   Tube-side Re = %.0f: %s flow → %s heat transfer\n', results.Re_t, ...
        iif(results.Re_t > 10000, 'Turbulent', iif(results.Re_t > 2300, 'Transitional', 'Laminar')), ...
        iif(results.Re_t > 10000, 'Excellent', iif(results.Re_t > 2300, 'Good', 'Poor')));
fprintf(fid, '   Shell-side Re = %.0f: %s flow → %s heat transfer\n', results.Re_s, ...
        iif(results.Re_s > 10000, 'Turbulent', iif(results.Re_s > 2300, 'Transitional', 'Laminar')), ...
        iif(results.Re_s > 10000, 'Excellent', iif(results.Re_s > 2300, 'Good', 'Poor')));
fprintf(fid, '\n');

fprintf(fid, '3. COMPACTNESS:\n');
volume = pi/4 * results.D_s^2 * params.L_tube;
fprintf(fid, '   Approximate volume: %.3f m³\n', volume);
fprintf(fid, '   Area density: %.1f m²/m³\n', results.A_o/volume);
fprintf(fid, '   Power density: %.1f kW/m³\n', params.Q/1000/volume);
fprintf(fid, '\n');

fprintf(fid, '4. OPTIMIZATION OPPORTUNITIES:\n');
fprintf(fid, '   - To reduce size: Increase U by improving h_i or h_s\n');
fprintf(fid, '   - To reduce tube-side ΔP: Increase D_internal or reduce N_p\n');
fprintf(fid, '   - To reduce shell-side ΔP: Increase baffle_spacing_ratio or p_f\n');
fprintf(fid, '   - To improve performance: Decrease p_f or increase l_f (watch ΔP!)\n');
fprintf(fid, '\n\n');

% Footer
fprintf(fid, '================================================================================\n');
fprintf(fid, '                         END OF DESIGN SPECIFICATION\n');
fprintf(fid, '================================================================================\n');

fclose(fid);

fprintf('Design specification document created: %s\n', filename);
fprintf('Document saved in: %s\n', pwd);
fprintf('\n');

% Helper function for inline if
function out = iif(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end
