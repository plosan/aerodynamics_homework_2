clear;
close all;
clc;

%% 1. DATA LOADING
% 1.1 MH 60 Data
data_mh60_re100k = load('data/mh60_re100k').data;
data_mh60_re200k = load('data/mh60_re200k').data;
data_mh60_re400k = load('data/mh60_re400k').data;

% 1.2 MH 61 Data
data_mh61_re100k = load('data/mh61_re100k').data;
data_mh61_re200k = load('data/mh61_re200k').data;
data_mh61_re400k = load('data/mh61_re400k').data;
data_mh61_re1000k = load('data/mh61_re1000k').data;


%% 2. MH 60 Lift curve

% 2.1 Zero lift angle and lift coefficient computation
p_mh60_re200k = polyfit(data_mh60_re200k(1:29,1), data_mh60_re200k(1:29,2), 1);
Cl_alpha_mh60_re200k = p_mh60_re200k(1);
alpha_l0_mh60_re200k = -p_mh60_re200k(2)/p_mh60_re200k(1);

p_mh60_re400k = polyfit(data_mh60_re400k(1:26,1), data_mh60_re400k(1:26,2), 1);
Cl_alpha_mh60_re400k = p_mh60_re400k(1);
alpha_l0_mh60_re400k = -p_mh60_re400k(2)/p_mh60_re400k(1);

fprintf("MH 60\nLift curve\n");
fprintf("Re = 200k%10s = %.4f%10s = %.2fº\n", ...
    "Cl_alpha", Cl_alpha_mh60_re200k, "alpha_l0", alpha_l0_mh60_re200k);
fprintf("Re = 400k%10s = %.4f%10s = %.2fº\n", ...
    "Cl_alpha", Cl_alpha_mh60_re400k, "alpha_l0", alpha_l0_mh60_re400k);

% 2.2 Legend
legend_str = cell(1,3);
legend_str(1) = {sprintf("$\\mathrm{Re} = 100k$")};
legend_str(2) = {sprintf("$\\mathrm{Re} = 200k$")};
legend_str(3) = {sprintf("$\\mathrm{Re} = 400k$")};

% 2.3 Plot
figure();
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{MH 60 Coeficiente de sustentaci\'on}");
scatter(data_mh60_re100k(:,1), data_mh60_re100k(:,2), 20, 'filled', 'r');
scatter(data_mh60_re200k(:,1), data_mh60_re200k(:,2), 20, 'filled', 'g');
scatter(data_mh60_re400k(:,1), data_mh60_re400k(:,2), 20, 'filled', 'b');
plot(data_mh60_re200k(1:29,1), polyval(p_mh60_re200k, data_mh60_re200k(1:29,1)), 'g');
plot(data_mh60_re400k(1:26,1), polyval(p_mh60_re400k, data_mh60_re400k(1:26,1)), 'b');
xlabel("\'Angulo de ataque $\left( \mathrm{deg} \right)$");
ylabel("Coeficiente de sustentaci\'on");
set(gcf, 'units', 'centimeters', 'position', [0,1,18,15]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
legend(legend_str, 'Location', 'Northwest');
grid on;
grid minor;
box on;
hold off;

%% 3. MH 60 Drag curve

% 3.1 Polar curve coefficients computation
p_mh60_re200k = polyfit(data_mh60_re200k(1:29,2).^2, data_mh60_re200k(1:29,3), 2);
p_mh60_re200k(2) = 0;
k_mh60_re200k = p_mh60_re200k(1);
Cd0_mh60_re200k = p_mh60_re200k(3);

p_mh60_re400k = polyfit(data_mh60_re400k(1:26,2).^2, data_mh60_re400k(1:26,3), 2);
p_mh60_re400k(2) = 0;
k_mh60_re400k = p_mh60_re400k(1);
Cd0_mh60_re400k = p_mh60_re400k(3);

fprintf("Drag curve\n");
fprintf("Re = 200k%10s = %.4f%10s = %.4f\n", ...
    "k", k_mh60_re200k, "CD0", Cd0_mh60_re200k);
fprintf("Re = 400k%10s = %.4f%10s = %.4f\n", ...
    "k", k_mh60_re400k, "CD0", Cd0_mh60_re400k);

% 3.2 Drag coefficient
figure();
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{MH 60 Coeficiente de resistencia aerodin\'amica}");
scatter(data_mh60_re100k(:,2), data_mh60_re100k(:,3), 20, 'filled', 'r');
scatter(data_mh60_re200k(:,2), data_mh60_re200k(:,3), 20, 'filled', 'g');
scatter(data_mh60_re400k(:,2), data_mh60_re400k(:,3), 20, 'filled', 'b');
plot(data_mh60_re200k(1:29,2), polyval(p_mh60_re200k, data_mh60_re200k(1:29,2)), 'g');
plot(data_mh60_re400k(1:26,2), polyval(p_mh60_re400k, data_mh60_re400k(1:26,2)), 'b');
xlabel("Coeficiente de sustentaci\'on");
ylabel("Coeficiente de resistencia aerodin\'amica");
% ylim([0 2.5]);
% yticks([0:0.25:2.5]);
set(gcf, 'units', 'centimeters', 'position', [0,1,18,15]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
legend(legend_str, 'Location', 'Northwest');
grid on;
grid minor;
box on;
hold off;

% % 2.4 Moment coefficient
% figure();
% hold on;
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% title("\textbf{MH 60 Coeficiente de momento libre}");
% plot(data_mh60_re100k(:,2), data_mh60_re100k(:,4), 'r');
% plot(data_mh60_re200k(:,2), data_mh60_re200k(:,4), 'g');
% plot(data_mh60_re400k(:,2), data_mh60_re400k(:,4), 'b');
% xlabel("Coeficiente de sustentaci\'on");
% ylabel("Coeficiente de momento libre");
% % ylim([0 2.5]);
% % yticks([0:0.25:2.5]);
% set(gcf, 'units', 'centimeters', 'position', [0,1,18,15]);
% % set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% % set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
% legend(legend_str, 'Location', 'Northwest');
% grid on;
% grid minor;
% box on;
% hold off;

%% 3. MH 61 PLOTS

legend_str(end+1) = {sprintf("$\\mathrm{Re} = 1000k$")};

% 3.1 Lift coefficient
figure();
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{MH 61 Coeficiente de sustentaci\'on}");
plot(data_mh61_re100k(:,1), data_mh61_re100k(:,2), 'r');
plot(data_mh61_re200k(:,1), data_mh61_re200k(:,2), 'g');
plot(data_mh61_re400k(:,1), data_mh61_re400k(:,2), 'b');
plot(data_mh61_re1000k(:,1), data_mh61_re1000k(:,2), 'm');
xlabel("\'Angulo de ataque $\left( \mathrm{deg} \right)$");
ylabel("Coeficiente de sustentaci\'on");
set(gcf, 'units', 'centimeters', 'position', [18,1,18,15]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
legend(legend_str, 'Location', 'Northwest');
grid on;
grid minor;
box on;
hold off;

% 3.2 Drag coefficient
figure();
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{MH 61 Coeficiente de resistencia aerodin\'amica}");
plot(data_mh61_re100k(:,2), data_mh61_re100k(:,3), 'r');
plot(data_mh61_re200k(:,2), data_mh61_re200k(:,3), 'g');
plot(data_mh61_re400k(:,2), data_mh61_re400k(:,3), 'b');
plot(data_mh61_re1000k(:,2), data_mh61_re1000k(:,3), 'm');
xlabel("Coeficiente de sustentaci\'on");
ylabel("Coeficiente de resistencia aerodin\'amica");
set(gcf, 'units', 'centimeters', 'position', [18,1,18,15]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
legend(legend_str, 'Location', 'Northwest');
grid on;
grid minor;
box on;
hold off;

% % 3.3 Moment coefficient
% figure();
% hold on;
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% title("\textbf{MH 61 Coeficiente de momento libre}");
% plot(data_mh61_re100k(:,2), data_mh61_re100k(:,4), 'r');
% plot(data_mh61_re200k(:,2), data_mh61_re200k(:,4), 'g');
% plot(data_mh61_re400k(:,2), data_mh61_re400k(:,4), 'b');
% plot(data_mh61_re1000k(:,2), data_mh61_re1000k(:,4), 'm');
% xlabel("Coeficiente de sustentaci\'on");
% ylabel("Coeficiente de momento libre");
% set(gcf, 'units', 'centimeters', 'position', [18,1,18,15]);
% % set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% % set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
% legend(legend_str, 'Location', 'Northwest');
% grid on;
% grid minor;
% box on;
% hold off;

%% 9. OTHER PLOTS
legend_str = cell(1,2);
legend_str(1) = {sprintf("MH 60")};
legend_str(2) = {sprintf("MH 61")};

mh60_ClCd_2 = data_mh60_re200k(:,2)./data_mh60_re200k(:,3);
mh61_ClCd_2 = data_mh61_re200k(:,2)./data_mh61_re200k(:,3);

mh60_ClCd_4 = data_mh60_re400k(:,2)./data_mh60_re400k(:,3);
mh61_ClCd_4 = data_mh61_re400k(:,2)./data_mh61_re400k(:,3);

figure();
hold on;
title("\textbf{Eficiencia aerodin\'amica para $\mathrm{Re} = 200k$}");
plot(data_mh60_re200k(:,1), mh60_ClCd_2, 'b');
plot(data_mh61_re200k(:,1), mh61_ClCd_2, 'r');
xlabel("\'Angulo de ataque $\left( \mathrm{deg} \right)$");
ylabel("Eficiencia aerodin\'amica");
set(gcf, 'units', 'centimeters', 'position', [0,1,18,15]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
legend(legend_str);
grid on;
grid minor;
box on;
hold off;

figure();
hold on;
title("\textbf{Eficiencia aerodin\'amica para $\mathrm{Re} = 400k$}");
plot(data_mh60_re400k(:,1), mh60_ClCd_4, 'b');
plot(data_mh61_re400k(:,1), mh61_ClCd_4, 'r');
xlabel("\'Angulo de ataque $\left( \mathrm{deg} \right)$");
ylabel("Eficiencia aerodin\'amica");
set(gcf, 'units', 'centimeters', 'position', [18,1,18,15]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
% set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
legend(legend_str);
grid on;
grid minor;
box on;
hold off;




