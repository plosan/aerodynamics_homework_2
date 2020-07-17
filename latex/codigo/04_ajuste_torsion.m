%% ETIP COMPUTATION
N = 100; % number of panels along the span
% 1 Constants
g0 = 9.80665;   % Gravity acceleration at SL        [m/s^2]
rho = 1.225;    % Air density at SL                 [kg/m^3]
V = 100/3.6;    % Flight speed                      [m/s]
% 2 Flying wing properties
WS = 20;                            % Wing load     [kg/m^2]
CL_design = g0*WS/(0.5*rho*V^2);    % Design CL     [1]
% 3 Solution vectors
CM_cg = zeros(1, length(ETIP));     % Moment coef. around CG            [1]
alpha_d = zeros(1, length(ETIP));   % Angle of attack in design cond    [deg]
% 4 Compute CG Moment coefficient and alpha design for several ETIP
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
end
% 5 Compute ETIP to trim the wing
i = 1;  found = 0;
ETIP_trim = 0;  CM_cg_trim = 0;     alpha_trim = 0;
while (i <= length(ETIP)-1) && (found == 0)
    if CM_cg(i)*CM_cg(i+1) < 0 % CM_cg = 0 for some ETIP between i and i+1
        x = CM_cg(i+1)/(CM_cg(i+1)-CM_cg(i));
        ETIP_trim = x*ETIP(i) + (1-x)*ETIP(i+1);
        CM_cg_trim = x*CM_cg(i) + (1-x)*CM_cg(i+1);
        alpha_trim = x*alpha_d(i) + (1-x)*alpha_d(i+1);
        found = 1;
    end
    i = i + 1;
end








