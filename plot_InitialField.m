plot_settings(1).variable = 'p';
plot_settings(2).variable = 'M';
plot_settings(3).variable = 's';
plot_settings(4).variable = 'sigma';
plot_settings(5).variable = 'theta';
plot_settings(6).variable = 'u';
plot_settings(1).level = linspace(0, 1, 21);
plot_settings(1).level_c = linspace(0, 1, 11);
plot_settings(2).level = linspace(60, 66, 31);
plot_settings(2).level_c = linspace(60, 66, 16);
plot_settings(3).level = linspace(-1, 1, 21);
plot_settings(3).level_c = linspace(-1, 1, 11);
plot_settings(5).level = 300:0.5:311;
plot_settings(5).level_c = 300:1:311;
plot_settings(6).level = -1:0.1:1;
plot_settings(6).level_c = -1:0.2:1;

CustomColour

figure(1)
clf
% Pressure
subplot(2, 2, 1)
hold on
contourf(phi, th_pad, Pres(:, :).', ...
    plot_settings(1).level, "LineStyle", "none")
contour(phi, th_pad, Pres(:, :).', ...
    plot_settings(1).level_c, "LineColor", "k")
hold off
xlim([phis, phin])
ylim([thetab, thetat])
colorbar
box on

% M
subplot(2, 2, 2)
hold on
contourf(phi(2:2:end - 1), th_pad, X(2:2:end - 1, :).', ...
    plot_settings(2).level, "LineStyle", "none")
contour(phi(2:2:end - 1), th_pad, X(2:2:end - 1, :).', ...
    plot_settings(2).level_c, "LineColor", "k")
hold off
xlim([phis, phin])
ylim([thetab, thetat])
colorbar

% s
subplot(2, 2, 3)
hold on
contourf(phi(1:2:end), th_pad, X(1:2:end, :).', ...
    plot_settings(3).level, "LineStyle", "none")
contour(phi(1:2:end), th_pad, X(1:2:end, :).', ...
    plot_settings(3).level_c, "LineColor", "k")
hold off
xlim([phis, phin])
ylim([thetab, thetat])
colormap(ccmap.seismic)
colorbar

% sigma
subplot(2, 2, 4)
contourf(phi, th_pad, Sigma(:, :).')
xlim([phis, phin])
ylim([thetab, thetat])
colorbar
