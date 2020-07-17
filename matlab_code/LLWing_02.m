% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodinàmica ESEIAAT-UPC
% e.ortega@upc.edu
% -------------------------------------------------------------------------     

clc; clear; close all;
format long;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1. INPUT DATA

% Airfoil data
load('data/mh60_re6.mat');
load('data/mh60_re10.mat');

% Lift coef. vs. alpha curve regression
poly_cl6 = polyfit(mh60_re6(1,1:10), mh60_re6(2,1:10), 1);
poly_cl10 = polyfit(mh60_re10(1,1:11), mh60_re10(2,1:11), 1);
alpha_l0_6 = roots(poly_cl6);
alpha_l0_10 = roots(poly_cl10);

% Drag coef. vs Lift coef curve regression
poly_p6 = polyfit(mh60_re6(2,2:11).^2, mh60_re6(3,2:11), 1);
poly_p10 = polyfit(mh60_re10(2,3:11).^2, mh60_re10(3,3:11), 1);

% Wing planform adimensional parameters (assumes planar wing)
AR = 21.3;   % aspect ratio
TR = 1/5.55;   % taper ratio    
DE25 = 16.50; % sweep angle at c/4 (deg)
ETIP = -7.1; % tip twist (deg, negative for washout)

% Wing planform dimensional parameters
x_cg = 1.38;    % CG location (behing root section LE)  [m]
cr = 1.55;      % Root chord                            [m]
b = 20;         % Wing span                             [m]

% Compute Stability margin (%)
sm = computeStabilityMargin(TR, DE25, x_cg, cr, b, 0);

% Plot Stability margin vs. Wing sweep
% DE25_vec = 12:1e-1:89;
% sm = zeros(1, length(DE25_vec));
% for i = 1:length(DE25_vec)
%     sm(i) = computeStabilityMargin(TR, DE25_vec(i), x_cg, cr, b, 0);
% end
% plotStabilityMargin(DE25_vec, sm, 1);

% Sections data (uses linear interpolation between root and tip)
% A0p = [alpha_l0_10 alpha_l0_6]; % root and tip section zero-lift angles (deg)
A0p = [-0.5 -0.5]; % root and tip section zero-lift angles (deg)
CM0p = [0.0140 0.0140]; % root and tip section free moments (deg)
CDP = [poly_p10(2) poly_p10(1);   % root section CD0 & k  (parabolic polar)
       poly_p6(2)  poly_p6(1)];  % tip section CD0 & k 

% Flap/aileron (symmetrical deflection)
YF_pos = [0.0 0.0];    % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.0;         % flap_chord/chord ratio
DE_flap = 0.0;          % flap deflection (deg, positive:down)
FlapCorr = 0.0;         % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)
N = 25; % number of panels along the span
ALPHA = [-2 10]; % angles of attack for analysis (deg) 0


%% 2. LIFTING LINE SOLUTION

% Wing discretization (lenghts are dimensionless with the wing span)
[c4nods, c75nods, chord, s_pan, h, Cm0_y, normals, mac, S] = ...
    geo(AR, TR, N, DE25, ETIP, A0p, CM0p, CDP, YF_pos, CF_ratio, ...
    DE_flap, FlapCorr); 

% Assembly of the influence coefficients matrix (needed once)
[inv_A, wake_len] = infcoeff(N, c4nods, c75nods, normals, h);

% Solve circulations for different angles of attack
[GAMMA, Ui, ncases] = getcirc(N, ALPHA, inv_A, normals);

% Loads calculation using plane Kutta-Joukowsky theorem (costlier, ...
% but general... and illustrative!)
[cl_local, force_coeff] = KuttaJoukowsky(N, c4nods, h, GAMMA, Ui, s_pan,...
    Cm0_y, chord, CDP, ncases, wake_len, S, mac, ALPHA);








