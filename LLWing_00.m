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
YF_pos = [0.0 0.0];     % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.0;         % flap_chord/chord ratio
DE_flap = 0.0;          % flap deflection (deg, positive:down)
FlapCorr = 0.0;         % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)
N = 1:1:200;            % number of panels along the span
ALPHA = [0 5];        % angles of attack for analysis (deg) 0

%% 2. CONVERGENCE ANALYSIS

CL_0 = zeros(1, length(N));
CM_LE_0 = zeros(1, length(N));
CL_5 = zeros(1, length(N));
CM_LE_5 = zeros(1, length(N));

for i = 1:length(N)    
    % Lifting Line Solution
    [c4nods, c75nods, chord, s_pan, h, Cm0_y, normals, mac, S] = ...
        geo(AR, TR, N(i), DE25, ETIP, A0p, CM0p, CDP, YF_pos, CF_ratio, ...
        DE_flap, FlapCorr); 
    [inv_A, wake_len] = infcoeff(N(i), c4nods, c75nods, normals, h);
    [GAMMA, Ui, ncases] = getcirc(N(i), ALPHA, inv_A, normals);
    [cl_local, force_coeff] = KuttaJoukowsky(N(i), c4nods, h, GAMMA, Ui, s_pan,...
    Cm0_y, chord, CDP, ncases, wake_len, S, mac, ALPHA);
    
    % CL, CM_LE at 0 and 5 deg
    CL_0(i) = force_coeff(7,1);
    CL_5(i) = force_coeff(7,2);
    CM_LE_0(i) = force_coeff(5,1);
    CM_LE_5(i) = force_coeff(5,2);
end

h = figure();
hold on;
title("\textbf{An\'alisis de convergencia $C_L$ y $C_{M,LE}$}");
plot(N, CL_0, 'r');
plot(N, CL_5, 'g');
plot(N, CM_LE_0, 'b');
plot(N, CM_LE_5, 'm');
xlabel("N\'umero de paneles");
legend("Coeficiente de sustentaci\'on, $\alpha = 0 ^\circ$", ...
    "Coeficiente de sustentaci\'on, $\alpha = 5 ^\circ$", ...
    "Coeficiente de momento LE, $\alpha = 0 ^\circ$", ...
    "Coeficiente de momento LE, $\alpha = 5 ^\circ$", 'NumColumns', 2);
set(gcf, 'units', 'centimeters', 'position', [0,1,18,11]);
grid on; grid minor; box on;    
hold off;

% Print
set(h, 'Units', 'Centimeters');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)])
print(h, 'plot/convergence_analysis', '-dpdf', '-r0');









