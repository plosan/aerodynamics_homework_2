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
ETIP = [-12:0.5:0]; % tip twist (deg, negative for washout)

% Wing planform dimensional parameters
x_cg = 1.38;    % CG location (behing root section LE)  [m]
cr = 1.55;      % Root chord                            [m]
b = 20;         % Wing span                             [m]

% Compute Stability margin (%)
[sm, x_ac, mean_ac] = computeStabilityMargin2(TR, DE25, x_cg, cr, b, 0);

% Plot Stability margin vs. Wing sweep
% DE25_vec = 12:1e-1:89;
% sm = zeros(1, length(DE25_vec));
% for i = 1:length(DE25_vec)
%     sm(i) = computeStabilityMargin(TR, DE25_vec(i), x_cg, cr, b, 0);
% end
% plotStabilityMargin(DE25_vec, sm, 1);

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
N = 100; % number of panels along the span


%% 2. ETIP COMPUTATION

% 2.1 Constants
g0 = 9.80665;   % Gravity acceleration at SL        [m/s^2]
rho = 1.225;    % Air density at SL                 [kg/m^3]
V = 100/3.6;    % Flight speed                      [m/s]

% 2.2 Flying wing properties
WS = 20;                            % Wing load     [kg/m^2]
CL_design = g0*WS/(0.5*rho*V^2);    % Design CL     [1]

% 2.3 Solution vectors
CM_cg = zeros(1, length(ETIP));     % Moment coef. around CG            [1]
alpha_d = zeros(1, length(ETIP));   % Angle of attack in design cond    [deg]

% 2.4 Compute CG Moment coefficient and alpha design for several ETIP
fprintf("%6s %12s %12s %12s %12s\n", "ETIP", "Cm_cg", "alpha_d", "x_ac", "x_cp");
for i = 1:length(ETIP)
    % Angles of attack for analysis
    ALPHA = -2:1:10;
    % Lifting line solution
    [c4nods, c75nods, chord, s_pan, h, Cm0_y, normals, mac, S] = ...
        geo(AR, TR, N, DE25, ETIP(i), A0p, CM0p, CDP, YF_pos, CF_ratio, ...
        DE_flap, FlapCorr); 
    [inv_A, wake_len] = infcoeff(N, c4nods, c75nods, normals, h);
    [GAMMA, Ui, ncases] = getcirc(N, ALPHA, inv_A, normals);
    [cl_local, force_coeff] = KuttaJoukowsky(N, c4nods, h, GAMMA, Ui, s_pan,...
        Cm0_y, chord, CDP, ncases, wake_len, S, mac, ALPHA);
    % Linear regression of CL and CM curves
    poly_CL = polyfit(ALPHA, force_coeff(7,:), 1);
    poly_CM = polyfit(force_coeff(7,:), force_coeff(5,:), 1);
    % Compute CM_LE and alpha in design condition
    CM_LE_design = polyval(poly_CM, CL_design);
    alpha_d(i) = (CL_design - poly_CL(2))/poly_CL(1);
    % CG Moment coefficient at design condition
    x_ac = -mac*b*poly_CM(1);                           % Wing aerodynamic center   [adim]
    x_cp = -mac*b*CM_LE_design/CL_design;               % Wing center of pressure   [adim]
    CM0 = CL_design*(x_ac - x_cp)/(mac*b);              % Wing free moment          [adim]
    CM_cg(i) = CM0 - CL_design*(x_ac - x_cg)/(mac*b);   % CG moment                 [adim]
    fprintf("%6.1f %12.4f %12.4f %12.4f %12.4f\n", ETIP(i), CM_cg(i), alpha_d(i), x_ac, x_cp);
end

% 2.5 Compute ETIP to trim the wing
i = 1;
found = 0;
ETIP_trim = 0;
CM_cg_trim = 0;
alpha_trim = 0;
while (i <= length(ETIP)-1) && (found == 0)
    if CM_cg(i)*CM_cg(i+1) < 0
        x = CM_cg(i+1)/(CM_cg(i+1)-CM_cg(i));
        ETIP_trim = x*ETIP(i) + (1-x)*ETIP(i+1);
        CM_cg_trim = x*CM_cg(i) + (1-x)*CM_cg(i+1);
        alpha_trim = x*alpha_d(i) + (1-x)*alpha_d(i+1);
        found = 1;
    end
    i = i + 1;
end

% 2.6 Print solution
fprintf("%15s = %.4f %s\n", "ETIP_t", ETIP_trim, "deg");
fprintf("%15s = %.4f %s\n", "CM_cg_trim", CM_cg_trim, "");
fprintf("%15s = %.4f %s\n", "alpha_trim", alpha_trim, "deg");


%% 3. PLOT
% 3.1 Plot CM_cg and alpha_design versus ETIP
h = figure(1);
hold on;
title("\textbf{Coeficiente de momento CG, \'angulo de ataque en condici\'on de dise\~no}");
yyaxis left;
plot(ETIP, CM_cg, 'b');
scatter(ETIP_trim, 0, 20, 'b', 'filled');
ylabel("Coeficiente de momento CG");
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
set(gca,'ycolor','k');
yyaxis right;
ylabel("\'Angulo de ataque $\left( \mathrm{deg} \right)$");
plot(ETIP, alpha_d, 'r');
scatter(ETIP_trim, alpha_trim, 20, 'r', 'filled');
plot([ETIP_trim ETIP_trim], [ceil(max(alpha_d)) floor(min(alpha_d))], '--k');
xlabel("Torsi\'on en punta de ala $\left( \mathrm{deg} \right)$");
xticks([min(ETIP):1:max(ETIP)]);
legend("Coeficiente de momento CG", "Coeficiente de momento CG condici\'on dise\~no", ...
    "\'Angulo de ataque", "\'Angulo de ataque condici\'on dise\~no");
set(gca,'ycolor','k');
set(gcf, 'units', 'centimeters', 'position', [0,1,18,11]);
grid on; grid minor; box on;    
hold off;

% 3.2 Save plot as pdf
set(h, 'Units', 'Centimeters');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(h, 'plot/trim_wing', '-dpdf', '-r0');









