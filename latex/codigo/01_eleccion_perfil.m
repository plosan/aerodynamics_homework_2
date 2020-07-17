%% 1. DATA
% Airfoil MH 60 
mh60_6 = load('mh60_re6.mat').mh60_re6;     % Re = 6e5
mh60_10 = load('mh60_re10.mat').mh60_re10;  % Re = 1e6
% Airfoil MH 61
mh61_6 = load('mh61_re6.mat').mh61_re6;     % Re = 6e5
mh61_10 = load('mh61_re10.mat').mh61_re10;  % Re = 1e6
%% 2. MH 60 LIFT REGRESSION
poly_cl6 = polyfit(mh60_6(1,1:10), mh60_6(2,1:10), 1);     % Re = 6e5     
poly_cl10 = polyfit(mh60_10(1,1:11), mh60_10(2,1:11), 1);  % Re = 1e6
alpha_l0_6 = roots(poly_cl6);       % Zero lift angle for Re = 6e5
alpha_l0_10 = roots(poly_cl10);     % Zero lift angle for Re = 1e6
%% 3. MH 60 POLAR REGRESSION
poly_p6 = polyfit(mh60_6(2,2:11).^2, mh60_6(3,2:11), 1);    % Re = 6e5
poly_p10 = polyfit(mh60_10(2,3:11).^2, mh60_10(3,3:11), 1); % Re = 1e6
%% 4. MH 60 AND 61 AERODYNAMIC EFFICIENCY COMPARISON
eff_mh60_6 = mh60_6(2,:)./mh60_6(3,:);      % Re = 6e5
eff_mh60_10 = mh60_10(2,:)./mh60_10(3,:);   % Re = 1e6
eff_mh61_6 = mh61_6(2,:)./mh61_6(3,:);      % Re = 6e5
eff_mh61_10 = mh61_10(2,:)./mh61_10(3,:);   % Re = 1e6


