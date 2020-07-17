function sm = computeStabilityMargin(TR, DE25, x_cg, cr, b, plot_fig)
% Previous computations
ct = cr/TR;         % Tip chord 
DE25 = DE25*pi/180; % Conversion to radians
% Origins 
O = [0, 0];                         % c/4 of root chord
O_prime = [b/2, -b*tan(DE25)/2];    % c/4 of tip chord
% Semi-wing corners
A = O + [0, cr/4];              % Root section LE
C = O + [0, -3*cr/4];           % Root section TE
B = O_prime + [0, ct/4];        % Tip section LE
D = O_prime + [0, -3*ct/4];     % Tip section TE
% Centroid computation points
E = A + [0, ct];    % Root section LE + tip chord
G = C + [0, -ct];   % Root section TE - tip chord
F = B + [0, cr];    % Tip section LE + root chord
H = D + [0, -cr];   % Tip section TE - root chord
Z = A + [0, -x_cg]; % Gravity center location
% Lines equations
syms x;
AB = @(x) A(2) + x*(B(2)-A(2))/(B(1)-A(1)); % Wing's LE
CD = @(x) C(2) + x*(D(2)-C(2))/(D(1)-C(1)); % Wing's TE
EH = @(x) E(2) + x*(H(2)-E(2))/(H(1)-E(1)); % Line joining E and H
GF = @(x) G(2) + x*(F(2)-G(2))/(F(1)-G(1)); % Line joining G and F
OO = @(x) O(2) + x*(O_prime(2)-O(2))/(O_prime(1)-O(1)); % c/4 line
% Intersection between EH and GF
sol_x = vpasolve(E(2)+x*(H(2)-E(2))/(H(1)-E(1))-(G(2) + x*(F(2)-G(2))/(F(1)-G(1))), x);
sol_y = OO(sol_x);
% Mean aerodynamic chord
p1 = AB(sol_x);
p2 = CD(sol_x);
mac = p1-p2;
% Aerodynamic center coordinates
ACC = [0, sol_y];
% Stability margin
sm = 100*norm(ACC-Z)/mac;
% Plot
if plot_fig == 1 % Plot the geometrical method ...