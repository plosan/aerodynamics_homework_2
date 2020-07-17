function plotStabilityMargin(DE25_vec, sm, printpdf)

h = figure();
hold on;
title("\textbf{Margen de estabilidad}");
plot(DE25_vec, sm, 'b');
xlabel("\'Angulo de flecha $\left( \mathrm{deg} \right)$");
ylabel("Margen de estabilidad $\left( \% \right)$");
set(gcf, 'units', 'centimeters', 'position', [0,1,18,11]);
% set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
grid on; grid minor; box on;
hold off;

if printpdf == 1
    set(h, 'Units', 'Centimeters');
    pos = get(h, 'Position');
    set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
        'PaperSize',[pos(3), pos(4)])
    print(h,'plot/stability_margin','-dpdf','-r0');
end


end