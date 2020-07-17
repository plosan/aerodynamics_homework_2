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
ALPHA = [0 10 20 30]; 

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

% 2.4 Numerical integral
y = linspace(-1/2, 1/2, N);
int_Cla = 0;
int_Clb = 0;
for i = 1:N-1
    % Average values
    avg_Cla = (Cla(i)*chord(i)*b + Cla(i+1)*chord(i+1)*b)/2;
    avg_Clb = b*(Clb(i)*chord(i)*b + Clb(i+1)*chord(i+1)*b)/2;
    int_Cla = int_Cla + avg_Cla*(b/N);
    int_Clb = int_Clb + avg_Clb*(b/N);
end
int_Cla = int_Cla/(S*b^2);
int_Clb = int_Clb/(S*b^2);
fprintf("%10s = %.9f\n", "int_Clb", int_Clb);
fprintf("%10s = %.9f\n", "int_Cla", int_Cla);

% 2.5 Plot
h = figure(1);
hold on;
title("\textbf{Distribuci\'on de sustentaci\'on b\'asica y adicional}");
plot(y, Clb, 'b');
plot(y, Cla, 'r');
xlabel("Envergadura adimensional");
ylabel("Coeficiente de sustentaci\'on");
legend("Sustentaci\'on b\'asica", "Sustentaci\'on adicional", ...
    'Location', 'East');
set(gcf, 'units', 'centimeters', 'position', [0,1,18,11]);
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.1f'));
grid on; grid minor; box on;    
hold off;

% 2.6 Save plot as pdf
set(h, 'Units', 'Centimeters');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
    'PaperSize',[pos(3), pos(4)]);
print(h, 'plot/lift_distribution', '-dpdf', '-r0');










