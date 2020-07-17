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
YF_pos = [0.0 0.0];    % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.0;         % flap_chord/chord ratio
DE_flap = 0.0;          % flap deflection (deg, positive:down)
FlapCorr = 0.0;         % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)
N = 200; % number of panels along the span


%% 2. COMPUTATION OF CL AND CL DISTRIBUTION FOR TWO ANGLES OF ATTACK

% 2.1 Angles of attack for analysis
ALPHA = [0 10]; 

% 2.2 Lifting line solution
[c4nods, c75nods, chord, s_pan, h, Cm0_y, normals, mac, S] = ...
    geo(AR, TR, N, DE25, ETIP, A0p, CM0p, CDP, YF_pos, CF_ratio, ...
    DE_flap, FlapCorr); 
[inv_A, wake_len] = infcoeff(N, c4nods, c75nods, normals, h);
[GAMMA, Ui, ncases] = getcirc(N, ALPHA, inv_A, normals);
[cl_local, force_coeff] = KuttaJoukowsky(N, c4nods, h, GAMMA, Ui, s_pan,...
    Cm0_y, chord, CDP, ncases, wake_len, S, mac, ALPHA);

% 2.3 Computation of lift coefficient distributions
Cla = (cl_local(:,1) - cl_local(:,2))/(force_coeff(7,1) - force_coeff(7,2));
Clb = (cl_local(:,2)*force_coeff(7,1) - cl_local(:,1)*force_coeff(7,2))/(force_coeff(7,1) - force_coeff(7,2));

% 2.4 CL vs. Alpha regression
poly_CL = polyfit(ALPHA, force_coeff(7,:), 1);

%% 3. COMPUTATION OF WING'S CL_MAX

% 3.1 Section and wing properties
rt_chord = [max(chord) min(chord)]; % Root and tip chords for discretization
Cl_max = [1.1481 1.0932];           % Root and tip sections max Cl

% 3.2 Linear interpolation
x = (chord-rt_chord(2))/(rt_chord(1)-rt_chord(2)); 
Cl_max_panel = x*Cl_max(1) + (1-x)*Cl_max(2);

% 3.3 Constants and properties
g0 = 9.80665;       % Gravity acceleration at SL    [m/s^2]
rho = 1.225;        % Air density at SL             [kg/m^3]
WS = 20;            % Wing load                     [kg/m^2]
W = WS*S*b^2*g0;    % Airplane weight               [kg]

% 3.4 Computation of stall CL for each wing section
CL_ys = (Cl_max_panel - Clb)./Cla;

% 3.5 Computation of stall CL, stall angle of attack and stall speed
CL_stall = min(CL_ys);                              % Stall CL              [1]
alpha_stall = (CL_stall - poly_CL(2))/poly_CL(1);   % Stall angle of attack [deg]
V_stall = sqrt(2*W/(rho*S*b^2*CL_stall));           % Stall speed           [m/s]

i_stall = 0;
found = 0;
while (i_stall <= N) && (found == 0)
    i_stall = i_stall + 1;
    if CL_ys(i_stall) == CL_stall
        found = 1;
    end
end

% 3.6 Print results
fprintf("%15s = %.4f %s\n", "CL_stall", CL_stall, "");
fprintf("%15s = %.4f %s\n", "alpha_stall", alpha_stall, "deg");
fprintf("%15s = %.4f %s\n", "V_stall", V_stall, "m/s");
fprintf("%15s = %.4f %s\n", "V_stall", V_stall*3.6, "km/h");

% 3.7 Plot stall Cl distribution
Cl_stall = Clb + Cla*CL_stall;
y = linspace(-1/2, 1/2, N);

scatter_x = [-1/2+i_stall/N, +1/2-i_stall/N];
scatter_y = [CL_stall CL_stall];

h = figure();
hold on;
title("\textbf{Distribuciones de coeficiente de sustentaci\'on en p\'erdida}");
plot(y, Cl_stall, 'b');
plot(y, CL_ys, 'r');
scatter(scatter_x, scatter_y, 'r', 'filled');
xlabel("Envergadura adimensional");
ylabel("Coeficiente de sustentaci\'on");
yticks([0.5:0.25:2]);
set(gcf, 'units', 'centimeters', 'position', [0,1,18,11]);
legend("Coeficiente de sustentaci\'on local en p\'erdida", ...
    "Coeficiente de sustentaci\'on de ala $\left( C_L \right)$ necesario para p\'erdida", ...
    "Secciones en p\'erdida");
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
grid on; grid minor; box on;    
hold off;

% 3.8 Save plot as pdf
set(h, 'Units', 'Centimeters');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(h, 'plot/stall_lift_distribution', '-dpdf', '-r0');





























