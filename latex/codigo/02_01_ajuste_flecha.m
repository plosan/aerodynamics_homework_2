%% COMPUTE STABILITY MARGIN
% Wing planform dimensional parameters
x_cg = 1.38;    % CG location (behing root section LE)  [m]
cr = 1.55;      % Root chord                            [m]
b = 20;         % Wing span                             [m]
% Function to compute the stability margin (%)
sm = computeStabilityMargin(TR, DE25, x_cg, cr, b, 0);

