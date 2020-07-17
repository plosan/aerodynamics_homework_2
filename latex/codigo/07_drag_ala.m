%% 2. COMPUTATION OF AERODYNAMIC COEFFICIENTS FOR SEVERAL ANGLES OF ATTACK
N = 200; % number of panels along the span
% 2.1 Angles of attack for analysis
ALPHA = -2:1:12; 
% 2.2 Lifting line solution
[c4nods, c75nods, chord, s_pan, h, Cm0_y, normals, mac, S] = ...
    geo(AR, TR, N, DE25, ETIP, A0p, CM0p, CDP, YF_pos, CF_ratio, ...
    DE_flap, FlapCorr); 
[inv_A, wake_len] = infcoeff(N, c4nods, c75nods, normals, h);
[GAMMA, Ui, ncases] = getcirc(N, ALPHA, inv_A, normals);
[cl_local, force_coeff] = KuttaJoukowsky(N, c4nods, h, GAMMA, Ui, s_pan,...
    Cm0_y, chord, CDP, ncases, wake_len, S, mac, ALPHA);
%% 3 CURVES REGRESSION
% 3.1 Data
CL = force_coeff(7,:);
CD = force_coeff(11,:);
% 3.2 Regressions
poly_CL = polyfit(ALPHA, CL, 1);    % CL vs alpha regression
alpha_L0 = roots(poly_CL);          % Zero-lift angle
poly_CD = polyfit(CL, CD, 2);       % CD vs CL regression











