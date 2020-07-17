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
mh60_re6 = load('data/mh60_re6.mat').mh60_re6;
mh60_re10 = load('data/mh60_re10.mat').mh60_re10;

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
ETIP = -4.3677; % tip twist (deg, negative for washout)

% Sections data (uses linear interpolation between root and tip)
A0p = [alpha_l0_10 alpha_l0_6]; % root and tip section zero-lift angles (deg)
CM0p = [0.0140 0.0140]; % root and tip section free moments (deg)
CDP = [poly_p10(2) poly_p10(1);   % root section CD0 & k  (parabolic polar)
       poly_p6(2)  poly_p6(1)];  % tip section CD0 & k 

% Flap/aileron (symmetrical deflection)
YF_pos = [0.0 0.0];    % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.0;         % flap_chord/chord ratio
DE_flap = 0.0;          % flap deflection (deg, positive:down)
FlapCorr = 0.0;         % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)
N = 200; % number of panels along the span


%% 2. COMPUTATION OF AERODYNAMIC COEFFICIENTS FOR SEVERAL ANGLES OF ATTACK

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
poly_CL = polyfit(ALPHA, CL, 1);
alpha_L0 = roots(poly_CL);

poly_CD = polyfit(CL, CD, 2);

% 3.3 Print results
fprintf("Lift curve\n");
fprintf("%15s = %.4f %s\n", "CL_alpha", poly_CL(1), "");
fprintf("%15s = %.4f %s\n", "CL0", poly_CL(2), "");
fprintf("%15s = %.4f %s\n", "alpha_L0", alpha_L0, "deg");

fprintf("\nDrag curve\n");
fprintf("%15s = %.3e %s\n", "k", poly_CD(1), "");
fprintf("%15s = %.3e %s\n", "CD0", poly_CD(2), "");

% 3.3 CL vs alpha plot
leg_CL = cell(1,2);
leg_CL(1) = {sprintf("%s", "Simulaci\'on num\'erica")};
leg_CL(2) = {sprintf("$C_L = %.4f \\left( \\alpha + %.4f \\right)$", poly_CL(1), abs(poly_CL(2)))};

alpha_plot = linspace(min(ALPHA), max(ALPHA), 100);
CL_plot = polyval(poly_CL, alpha_plot);

figure(1);
hold on;
title("\textbf{Coeficiente de sustentaci\'on}");
plot(alpha_plot, CL_plot, 'b');
scatter(ALPHA, CL, 20, 'r', 'filled');
xlabel("\'Angulo de ataque");
ylabel("Coeficiente de sustentaci\'on");
set(gcf, 'units', 'centimeters', 'position', [0,1,18,11]);
legend(leg_CL, 'Location', 'Northwest');
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
grid on; grid minor; box on;    
hold off;

% 3.4 CD vs CL plot
leg_CD = cell(1,2);
leg_CD(1) = leg_CL(1);
leg_CD(2) = {sprintf("$C_D = 10^{-4} \\left( %.1f \\, {C_L}^2 + %.2f \\right)$", 1e4*poly_CD(1), 1e4*poly_CD(3))};
for i = 1:length(poly_CD)
    fprintf("%3d\t%.3f\n", i, 1e4*poly_CD(i));
end

poly_CD_plot = [poly_CD(1) 0 poly_CD(3)];

h = figure(2);
hold on;
title("\textbf{Resistencia aerodin\'amica}");
scatter(CL, CD, 20, 'r', 'filled');
plot(CL_plot, polyval(poly_CD_plot, CL_plot), 'b');
xlabel("Coeficiente de sustentaci\'on");
ylabel("Coeficiente de resistencia aerodin\'amica");
set(gcf, 'units', 'centimeters', 'position', [18,1,18,11]);
legend(leg_CD, 'Location', 'Northwest');
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.3f'));
grid on; grid minor; box on;    
hold off;

% 3.5 Save plot as pdf
set(h, 'Units', 'Centimeters');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(h, 'plot/drag_regression', '-dpdf', '-r0');











