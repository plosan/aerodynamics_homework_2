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

% Wing planform dimensional parameters
x_cg = 1.38;    % CG location (behing root section LE)  [m]
cr = 1.55;      % Root chord                            [m]
b = 20;         % Wing span                             [m]

% Sections data (uses linear interpolation between root and tip)
A0p = [alpha_l0_10 alpha_l0_6]; % root and tip section zero-lift angles (deg)
CM0p = [0.0140 0.0140]; % root and tip section free moments (deg)
CDP = [poly_p10(2) poly_p10(1);   % root section CD0 & k  (parabolic polar)
       poly_p6(2)  poly_p6(1)];  % tip section CD0 & k 

% Flap/aileron (symmetrical deflection)
YF_pos = 2*[3.75 7.75]/b;       % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2522/1.090;         % flap_chord/chord ratio
DE_flap = -20:5:20;                  % flap deflection (deg, positive:down)
FlapCorr = 0.8;                 % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)
N = 200; % number of panels along the span


%% 2. COMPUTATION OF CM_CG FOR SEVERAL FLAP DEFLECTIONS

% 2.1 Angles of attack for analysis
ALPHA = [-2 4 10]; 
legend_str = cell(length(DE_flap), 1);

% 2.2 Solution vectors
CL = zeros(length(DE_flap), length(ALPHA));
poly_CM = zeros(length(DE_flap), 2);
CM_cg = zeros(length(DE_flap), length(ALPHA));
poly_CM_cg = zeros(length(DE_flap), 2);

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
    CM_LE = force_coeff(5,:);
    CL(i,:) = force_coeff(7,:);
    poly_CM(i,:) = polyfit(CL(i,:), CM_LE, 1);
    x_ac = -mac*b*poly_CM(i,1);
    x_cp = -mac*b*CM_LE./CL(i,:);
    CM0 = CL(i,:).*(x_ac-x_cp)/(mac*b);
    CM_cg(i,:) = CM0 - CL(i,:)*(x_ac-x_cg)/(mac*b);
    poly_CM_cg(i,:) = polyfit(CL(i,:), CM_cg(i,:), 1);
    % Fill legend
    if DE_flap(i) > 0
        legend_str(i) = {sprintf("$\\delta_e = +%d ^\\circ$", DE_flap(i))};
    else
        legend_str(i) = {sprintf("$\\delta_e = %d ^\\circ$", DE_flap(i))};
    end
end

% 2.3 Load CD data
poly_CD = load("data/poly_CD.mat").poly_CD_plot;

CL_max_range = sqrt(poly_CD(3)/poly_CD(1));
CL_min_speed = sqrt(3)*CL_max_range;

fprintf("%15s = %.4f\n", "CL_max_range", CL_max_range);
fprintf("%15s = %.4f\n", "CL_min_speed", CL_min_speed);

fprintf("\n%15s = %.3e\n", "k", poly_CD(1));
fprintf("%15s = %.3e\n", "CD0", poly_CD(3));

% 2.3 Plot
h = figure(1);
hold on;
title("\textbf{Coeficiente de momento CG con deflexi\'on de flap}");
CL_plot = linspace(min(CL,[],'all'), max(CL,[],'all'), 10);

plot(CL_plot, polyval(poly_CM_cg(1,:), CL_plot), 'r');
plot(CL_plot, polyval(poly_CM_cg(2,:), CL_plot), 'g');
plot(CL_plot, polyval(poly_CM_cg(3,:), CL_plot), 'b');
plot(CL_plot, polyval(poly_CM_cg(4,:), CL_plot), 'm');
plot(CL_plot, polyval(poly_CM_cg(5,:), CL_plot), 'k');
plot(CL_plot, polyval(poly_CM_cg(6,:), CL_plot), '--r');
plot(CL_plot, polyval(poly_CM_cg(7,:), CL_plot), '--g');
plot(CL_plot, polyval(poly_CM_cg(8,:), CL_plot), '--b');
plot(CL_plot, polyval(poly_CM_cg(9,:), CL_plot), '--m');

scatter(CL(1,:), CM_cg(1,:), 20, 'r', 'filled');
scatter(CL(2,:), CM_cg(2,:), 20, 'g', 'filled');
scatter(CL(3,:), CM_cg(3,:), 20, 'b', 'filled');
scatter(CL(4,:), CM_cg(4,:), 20, 'm', 'filled');
scatter(CL(5,:), CM_cg(5,:), 20, 'k', 'filled');
scatter(CL(6,:), CM_cg(6,:), 20, 'r', 'filled');
scatter(CL(7,:), CM_cg(7,:), 20, 'g', 'filled');
scatter(CL(8,:), CM_cg(8,:), 20, 'b', 'filled');
scatter(CL(9,:), CM_cg(9,:), 20, 'm', 'filled');

ylim([-0.4 0.4]);

plot([CL_max_range CL_max_range], [-0.4 0.4], '--k');
plot([CL_min_speed CL_min_speed], [-0.4 0.4], '--k');

h1 = text(CL_max_range-0.05, 0.18, "M\'ax. alcance");
set(h1, 'Rotation', 90);

h2 = text(CL_min_speed-0.05, 0.18, "M\'in. velocidad");
set(h2, 'Rotation', 90);

xlabel("Coeficiente de sustentaci\'on");
ylabel("Coeficiente de momento CG");
legend(legend_str, 'NumColumns', 1, 'Location', 'Southwest');
set(gcf, 'units', 'centimeters', 'position', [0,1,18,11]);
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
grid on; grid minor; box on;    
hold off;

% 2.4 Save plot as pdf
set(h, 'Units', 'Centimeters');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(h, 'plot/cmcg_flap', '-dpdf', '-r0');











