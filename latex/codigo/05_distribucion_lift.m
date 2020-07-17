%% COMPUTATION OF CL AND CL DISTRIBUTION FOR TWO ANGLES OF ATTACK
% 1. Angles of attack for analysis
ALPHA = [0 10]; 
% 2. Lifting line solution
[c4nods, c75nods, chord, s_pan, h, Cm0_y, normals, mac, S] = ...
    geo(AR, TR, N, DE25, ETIP, A0p, CM0p, CDP, YF_pos, CF_ratio, ...
    DE_flap, FlapCorr); 
[inv_A, wake_len] = infcoeff(N, c4nods, c75nods, normals, h);
[GAMMA, Ui, ncases] = getcirc(N, ALPHA, inv_A, normals);
[cl_local, force_coeff] = KuttaJoukowsky(N, c4nods, h, GAMMA, Ui, s_pan,...
    Cm0_y, chord, CDP, ncases, wake_len, S, mac, ALPHA);
% 2. Computation of lift coefficient distributions
Cla = (cl_local(:,1) - cl_local(:,2))/(force_coeff(7,1) - force_coeff(7,2));
Clb = (cl_local(:,2)*force_coeff(7,1) - cl_local(:,1)*force_coeff(7,2))/(force_coeff(7,1) - force_coeff(7,2));
% 3. Plot

























