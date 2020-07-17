%% COMPUTATION OF WING'S CL_MAX
% 1. Section and wing properties, some obtained from LL solution in the
% computation of Cl distributions
rt_chord = [max(chord) min(chord)]; % Root and tip chords for discretization
Cl_max = [1.1481 1.0932];           % Root and tip sections max Cl
% 2. Linear interpolation
x = (chord-rt_chord(2))/(rt_chord(1)-rt_chord(2));  % Interpolation variable
Cl_max_panel = x*Cl_max(1) + (1-x)*Cl_max(2);       % Cl_max for each panel
% 3. Constants and properties
g0 = 9.80665;       % Gravity acceleration at SL    [m/s^2]
rho = 1.225;        % Air density at SL             [kg/m^3]
WS = 20;            % Wing load                     [kg/m^2]
W = WS*S*b^2*g0;    % Airplane weight               [kg]
% 4. Computation of stall CL for each wing section
CL_ys = (Cl_max_panel - Clb)./Cla;
% 5. Computation of stall CL, stall angle of attack and stall speed
CL_stall = min(CL_ys);                              % Stall CL              [1]
alpha_stall = (CL_stall - poly_CL(2))/poly_CL(1);   % Stall angle of attack [deg]
V_stall = sqrt(2*W/(rho*S*b^2*CL_stall));           % Stall speed           [m/s]
% 6. Search the first panel in stall
i_stall = 0;    found = 0;
while (i_stall <= N) && (found == 0)
    i_stall = i_stall + 1;
    if CL_ys(i_stall) == CL_stall
        found = 1;
    end
end