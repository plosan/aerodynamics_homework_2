%% 1. INPUT DATA
% Wing planform adimensional parameters (assumes planar wing)
AR = 21.3;   % aspect ratio
TR = 1/5.55;   % taper ratio    
DE25 = 16.50; % sweep angle at c/4 (deg)
ETIP = -7.1; % tip twist (deg, negative for washout)
% Sections data (uses linear interpolation between root and tip)
A0p = [alpha_l0_10 alpha_l0_6]; % root and tip section zero-lift angles (deg)
CM0p = [0.0140 0.0140]; % root and tip section free moments (deg)
CDP = [poly_p10(2) poly_p10(1);   % root section CD0 & k  (parabolic polar)
       poly_p6(2)  poly_p6(1)];  % tip section CD0 & k 
% Flap/aileron (symmetrical deflection)
YF_pos = [0.0 0.0];     % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.0;         % flap_chord/chord ratio
DE_flap = 0.0;          % flap deflection (deg, positive:down)
FlapCorr = 0.0;         % flap effectiviness (<=1)
% Simulation data (by the time being only longitudinal analysis)
N = 1:1:200;            % number of panels along the span
ALPHA = [0 5];        % angles of attack for analysis (deg) 0
%% 2. CONVERGENCE ANALYSIS
% 2.1 Solution matrices
CL_0 = zeros(1, length(N));     % Lift coefficient at 0 deg
CM_LE_0 = zeros(1, length(N));  % LE moment coefficient at 0 deg
CL_5 = zeros(1, length(N));     % Lift coefficient at 5 deg
CM_LE_5 = zeros(1, length(N));  % LE moment coefficient at 5 deg
% 2.2 Compute for every N
for i = 1:length(N)    
    % Lifting Line Solution
    [c4nods, c75nods, chord, s_pan, h, Cm0_y, normals, mac, S] = ...
        geo(AR, TR, N(i), DE25, ETIP, A0p, CM0p, CDP, YF_pos, CF_ratio, ...
        DE_flap, FlapCorr); 
    [inv_A, wake_len] = infcoeff(N(i), c4nods, c75nods, normals, h);
    [GAMMA, Ui, ncases] = getcirc(N(i), ALPHA, inv_A, normals);
    [cl_local, force_coeff] = KuttaJoukowsky(N(i), c4nods, h, GAMMA, Ui, s_pan,...
    Cm0_y, chord, CDP, ncases, wake_len, S, mac, ALPHA);    
    % Save computations
    CL_0(i) = force_coeff(7,1);
    CL_5(i) = force_coeff(7,2);
    CM_LE_0(i) = force_coeff(5,1);
    CM_LE_5(i) = force_coeff(5,2);
end
% 2.3 Plot









