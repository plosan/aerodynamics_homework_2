%% 1. INPUT DATA
% 1.1 Airfoil data
mh60_re6 = load('data/mh60_re6.mat').mh60_re6;
mh60_re10 = load('data/mh60_re10.mat').mh60_re10;
% 1.2 Lift coef. vs. alpha curve regression
poly_cl6 = polyfit(mh60_re6(1,1:10), mh60_re6(2,1:10), 1);
poly_cl10 = polyfit(mh60_re10(1,1:11), mh60_re10(2,1:11), 1);
alpha_l0_6 = roots(poly_cl6);
alpha_l0_10 = roots(poly_cl10);
% 1.3 Drag coef. vs Lift coef curve regression
poly_p6 = polyfit(mh60_re6(2,2:11).^2, mh60_re6(3,2:11), 1);
poly_p10 = polyfit(mh60_re10(2,3:11).^2, mh60_re10(3,3:11), 1);
% 1.4 Wing planform adimensional parameters (assumes planar wing)
AR = 21.3;   % aspect ratio
TR = 1/5.55;   % taper ratio    
DE25 = 16.50; % sweep angle at c/4 (deg)
ETIP = -4.3677; % tip twist (deg, negative for washout)
% 1.5 Sections data (uses linear interpolation between root and tip)
A0p = [alpha_l0_10 alpha_l0_6]; % root and tip section zero-lift angles (deg)
CM0p = [0.0140 0.0140]; % root and tip section free moments (deg)
CDP = [poly_p10(2) poly_p10(1);   % root section CD0 & k  (parabolic polar)
       poly_p6(2)  poly_p6(1)];  % tip section CD0 & k 
% 1.6 Flap/aileron (symmetrical deflection)
YF_pos = 2*[3.75 7.75]/b;       % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2522/1.090;         % flap_chord/chord ratio
DE_flap = -20:5:20;                  % flap deflection (deg, positive:down)
FlapCorr = 0.8;                 % flap effectiviness (<=1)

%% 2. COMPUTATION OF CM_CG FOR SEVERAL FLAP DEFLECTIONS
N = 200; % number of panels along the span
% 2.1 Angles of attack for analysis
ALPHA = [-2 4 10]; 
legend_str = cell(length(DE_flap), 1);  % Legend for plot
% 2.2 Solution vectors
CL = zeros(length(DE_flap), length(ALPHA));
poly_CM = zeros(length(DE_flap), 2);
CM_cg = zeros(length(DE_flap), length(ALPHA));
poly_CM_cg = zeros(length(DE_flap), 2);
% 2.3 Solve for each flap deflection angle
for i = 1:length(DE_flap)
    % Lifting line solution
    [c4nods, c75nods, chord, s_pan, h, Cm0_y, normals, mac, S] = ...
        geo(AR, TR, N, DE25, ETIP, A0p, CM0p, CDP, YF_pos, CF_ratio, ...
        DE_flap(i), FlapCorr); 
    [inv_A, wake_len] = infcoeff(N, c4nods, c75nods, normals, h);
    [GAMMA, Ui, ncases] = getcirc(N, ALPHA, inv_A, normals);
    [cl_local, force_coeff] = KuttaJoukowsky(N, c4nods, h, GAMMA, Ui, s_pan,...
        Cm0_y, chord, CDP, ncases, wake_len, S, mac, ALPHA);
    % Compute CM_cg
    CM_LE = force_coeff(5,:);   % Moment coefficient respecto to LE
    CL(i,:) = force_coeff(7,:); % Lift coefficient
    poly_CM(i,:) = polyfit(CL(i,:), CM_LE, 1);  % CM_LE vs CL regression
    x_ac = -mac*b*poly_CM(i,1);     % Aerodynamic center position
    x_cp = -mac*b*CM_LE./CL(i,:);   % Pressure center position
    CM0 = CL(i,:).*(x_ac-x_cp)/(mac*b); % Free moment coefficient
    CM_cg(i,:) = CM0 - CL(i,:)*(x_ac-x_cg)/(mac*b); % CG moment coefficient
    poly_CM_cg(i,:) = polyfit(CL(i,:), CM_cg(i,:), 1);  % CM_cg vs CL regression
    % Fill legend
    if DE_flap(i) > 0
        legend_str(i) = {sprintf("$\\delta_e = +%d ^\\circ$", DE_flap(i))};
    else
        legend_str(i) = {sprintf("$\\delta_e = %d ^\\circ$", DE_flap(i))};
    end
end
% 2.3 Load CD data and compute CL for max range and min speed conditions
poly_CD = load("data/poly_CD.mat").poly_CD_plot;
CL_max_range = sqrt(poly_CD(3)/poly_CD(1));
CL_min_speed = sqrt(3)*CL_max_range;









