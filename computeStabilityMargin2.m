function [sm, x_ac, mean_ac] = computeStabilityMargin2(TR, DE25, x_cg, cr, b, plot_fig)

ct = cr*TR;
DE25 = DE25*pi/180;

% Origins 
O = [0, 0];                         % c/4 of root chord
O_prime = [b/2, -b*tan(DE25)/2];    % c/4 of tip chord

% Semi-wing corners
A = O + [0, cr/4];
C = O + [0, -3*cr/4];
B = O_prime + [0, ct/4];
D = O_prime + [0, -3*ct/4];

% Centroid computation points
E = A + [0, ct];
G = C + [0, -ct];
F = B + [0, cr];
H = D + [0, -cr];

% Gravity center
Z = A + [0, -x_cg];

% Lines equations
syms x;
AB = @(x) A(2) + x*(B(2)-A(2))/(B(1)-A(1));
CD = @(x) C(2) + x*(D(2)-C(2))/(D(1)-C(1));
% EH = @(x) E(2) + x*(H(2)-E(2))/(H(1)-E(1));
% GF = @(x) G(2) + x*(F(2)-G(2))/(F(1)-G(1));
OO = @(x) O(2) + x*(O_prime(2)-O(2))/(O_prime(1)-O(1));

% Intersection between EH and GF
% int1 = @(x) E(2) + x*(H(2)-E(2))/(H(1)-E(1)) - (G(2) + x*(F(2)-G(2))/(F(1)-G(1)));
sol_x = vpasolve(E(2) + x*(H(2)-E(2))/(H(1)-E(1)) - (G(2) + x*(F(2)-G(2))/(F(1)-G(1))), x);
sol_y = OO(sol_x);
x_ac = cr/4 + abs(sol_y);

% Mean aerodynamic chord
p1 = AB(sol_x);
p2 = CD(sol_x);
mean_ac = p1-p2;

% Aerodynamic center coordinates
ACC = [0, sol_y];

% Stability margin
sm = 100*norm(ACC-Z)/mean_ac;

if plot_fig == 1
    figure();
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    hold on;
    axis equal;
    % Legend plots
    scatter(Z(1), Z(2), 50, 'g', 'filled');
    scatter(ACC(1), ACC(2), 50, 'm', 'filled');
    plot([sol_x sol_x], [p1 p2], 'c', 'LineWidth', 1);
    plot([O(1) O_prime(1)], [O(2) O_prime(2)], '--k', 'LineWidth', 1);
    % Wing plots
    plot([A(1) B(1)], [A(2) B(2)], 'k', 'LineWidth', 1);
    plot([B(1) D(1)], [B(2) D(2)], 'k', 'LineWidth', 1);
    plot([D(1) C(1)], [D(2) C(2)], 'k', 'LineWidth', 1);
    plot([A(1) C(1)], [A(2) C(2)], 'k', 'LineWidth', 1);
    % Auxiliary lines plots
    plot([E(1) H(1)], [E(2) H(2)], '--b');
    plot([G(1) F(1)], [G(2) F(2)], '--b');
    plot([A(1) E(1)], [A(2) E(2)], '--r');
    plot([C(1) G(1)], [C(2) G(2)], '--r');
    plot([B(1) F(1)], [B(2) F(2)], '--r');
    plot([D(1) H(1)], [D(2) H(2)], '--r');
    % Text plots
    text(A(1)-0.25, A(2), "$A$");
    text(B(1)+0.25, B(2), "$B$");
    text(C(1)-0.25, C(2), "$C$");
    text(D(1)+0.25, D(2), "$D$");
    text(E(1)-0.25, E(2), "$E$");
    text(F(1)+0.25, F(2), "$F$");
    text(G(1)-0.25, G(2), "$G$");
    text(H(1)+0.25, H(2), "$H$");
    xlim([-1 b/2+1]);
    ylim([-5 2]);
    xlabel("$y$");
    ylabel("$x$");
    legend("Gravity center", "Aerodynamic center", ...
        "Mean aerodynamic chord", "$c/4$ line");
    set(gcf, 'units', 'centimeters', 'position', [20,1,25,18]);
    grid on; grid minor; box on;    
    hold off;
end

end

