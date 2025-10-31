function results = calculateHeatExchanger(params)
% calculateHeatExchanger - Performs heat exchanger design calculation
%
% Input: params - structure containing all design parameters
% Output: results - structure containing calculated results
%
% Returns: U_o_calc, A_o, N_tt, DP_t_total, DP_s, converged, iterations

% ==========================================
% ========== EXTRACT PARAMETERS ============
% ==========================================

% Specifications
Q = params.Q;
m_hot = params.m_hot;
m_cold = params.m_cold;
T_hot1 = params.T_hot1;
T_cold1 = params.T_cold1;
F = params.F;

% Geometry
L_tube = params.L_tube;
D_internal = params.D_internal;
D_external = params.D_external;
k_tube = params.k_tube;

% Fins
l_f = params.l_f;
t_f = params.t_f;
p_f = params.p_f;
k_fin = params.k_fin;

% Layout
L_tp = params.L_tp;
N_p = params.N_p;
psi_n = params.psi_n;
C_1 = params.C_1;
L_bb = params.L_bb;
baffle_spacing_ratio = params.baffle_spacing_ratio;

% Fouling
Rf_i = params.Rf_i;
Rf_o = params.Rf_o;

% Calculation Parameters
U_initial_guess = params.U_initial_guess;
tolerance = params.tolerance;
j_f = params.j_f;

% Fluid Properties - Tube Side (Water)
rho_t = params.rho_t;
mu_t = params.mu_t;
k_t = params.k_t;
Cp_t = params.Cp_t;
Pr_t = (Cp_t * mu_t) / k_t;

% Fluid Properties - Shell Side (Air)
rho_s = params.rho_s;
mu_s = params.mu_s;
k_s = params.k_s;
Cp_s = params.Cp_s;
Pr_s = (Cp_s * mu_s) / k_s;

% ==========================================
% ===== ITERATIVE LOOP (THERMAL DESIGN) ====
% ==========================================

U_o_ass = U_initial_guess;  % Initial guess
max_iterations = 20;
converged = false;

for i = 1:max_iterations
    % --- Calculate Required Area
    T_hot2 = T_hot1 - Q / (m_hot * Cp_s);
    T_cold2 = T_cold1 + Q / (m_cold * Cp_t);
    dT1 = T_hot1 - T_cold2;
    dT2 = T_hot2 - T_cold1;
    DeltaT_LM = (dT1 - dT2) / log(dT1 / dT2);
    A_req = Q / (U_o_ass * DeltaT_LM * F);

    % --- Calculate Geometry
    r_f1 = D_external / 2;
    l_fc = l_f + t_f / 2;
    r_f2c = r_f1 + l_fc;
    A_f_one_fin = 2 * pi * (r_f2c^2 - r_f1^2);
    N_fin = L_tube / p_f;
    A_b = pi * D_external * (L_tube - N_fin * t_f);
    A_o_one_tube = (N_fin * A_f_one_fin) + A_b;

    N_tt = A_req / A_o_one_tube;

    % --- Calculate Tube-Side hi
    N_per_pass = N_tt / N_p;
    A_t_internal_one = (pi/4) * D_internal^2;
    A_t_flow = N_per_pass * A_t_internal_one;
    v_t = m_cold / (rho_t * A_t_flow);
    Re_t = (rho_t * v_t * D_internal) / mu_t;

    % Use correct exponent for HEATING (0.4 for heating)
    Nu_t = 0.023 * Re_t^0.8 * Pr_t^0.4;
    h_i = (Nu_t * k_t) / D_internal;

    % --- Calculate Shell-Side hs
    D_ctl = L_tp * (4*C_1*N_tt/(pi*(1-psi_n)))^(1/2);
    D_s = D_ctl + L_bb + D_external;
    L_B = baffle_spacing_ratio * D_s;
    A_face = L_tube * (D_ctl + D_external);
    Q_air = m_hot / rho_s;
    u_f = Q_air / A_face;
    u_max = u_f * L_tp / (L_tp - D_external);
    Re_s = (rho_s * u_max * D_external) / mu_s;

    % Nu_s correlation
    Nu_s = 0.134 * Re_s^0.681 * Pr_s^0.33 * ...
           ((p_f - t_f)/l_f)^0.2 * (p_f/t_f)^0.1134;
    h_s = (Nu_s * k_s) / D_external;

    % --- Calculate U_o_calc
    A_o = N_tt * A_o_one_tube;
    A_i = pi * D_internal * L_tube * N_tt;

    % Fin efficiency
    m = sqrt((2 * h_s) / (k_fin * t_f));
    eta_f = tanh(m * l_fc) / (m * l_fc);
    eta_o = 1 - (N_fin * A_f_one_fin / A_o) * (1 - eta_f);

    % Resistances (all referred to A_o)
    R_shell_conv = 1 / (eta_o * h_s);
    R_tube_conv = 1 / (h_i * (A_i / A_o));
    R_wall_cond = (D_external * log(D_external/D_internal)) / (2 * k_tube);
    R_fouling = (Rf_i * (A_o / A_i)) + (Rf_o / eta_o);

    R_total = R_shell_conv + R_tube_conv + R_wall_cond + R_fouling;
    U_o_calc = 1 / R_total;

    % --- Check Convergence
    error = abs(U_o_calc - U_o_ass);

    if error < tolerance
        converged = true;
        break;
    end

    % --- Update for Next Loop
    U_o_ass = U_o_calc;
end

% ============================================
% ==== PRESSURE DROP (HYDRAULIC DESIGN) =====
% ============================================

% --- Tube-Side DP
L_total = L_tube * N_p;
f_t = 0.0035 + 0.264 / Re_t^0.42;
DP_t_friction = 4 * f_t * (L_total / D_internal) * (rho_t * v_t^2 / 2);

K = 1.8 * N_p;
DP_t_return = K * (rho_t * v_t^2) / 2;
DP_t_total = DP_t_friction + DP_t_return;

% --- Shell-Side DP (Kern's Method)
A_s = ((L_tp - D_external) * D_s * L_B) / L_tp;
G_s = m_hot / A_s;
u_s = G_s / rho_s;
d_e = (1.10 / D_external) * (L_tp^2 - 0.917 * D_external^2);
Re_s_kern = (G_s * d_e) / mu_s;

% Use j_f from parameters (from chart lookup)
DP_s = 8 * j_f * (D_s / d_e) * (L_tube / L_B) * (rho_s * u_s^2 / 2);

% ==========================================
% ============ RETURN RESULTS ==============
% ==========================================

results.U_o_calc = U_o_calc;        % W/m2K
results.A_o = A_o;                  % m2
results.N_tt = N_tt;                % Number of tubes
results.DP_t_total = DP_t_total;    % Pa
results.DP_s = DP_s;                % Pa
results.converged = converged;      % Boolean
results.iterations = i;             % Number of iterations
results.D_s = D_s;                  % Shell diameter (m)
results.L_B = L_B;                  % Baffle spacing (m)
results.Re_t = Re_t;                % Tube-side Reynolds number
results.Re_s = Re_s;                % Shell-side Reynolds number

end
