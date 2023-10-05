export = false;

%% Intepolate to full domain
M_int = 2:2:Jt - 1;
s_int = 3:2:Jt - 2;

s_intp = X;
s_intp(M_int, :) = (X(M_int - 1, :) + X(M_int + 1, :)) / 2;
phi_intp = rad2deg(asin(s_intp));
Phi_intp = rad2deg(asin(YY));
M_intp = X;
M_intp(s_int, :) = (X(s_int - 1, :) + X(s_int + 1, :)) / 2;
Pres_intp = Pres;
Pres_intp(s_int, :) = (Pres(s_int - 1, :) + Pres(s_int + 1, :)) / 2;
Pres_intp = Pres_intp * pb / 100;
Pres0_intp = Pres0 * pb / 100;
u = Omega * Ae * ((sin(deg2rad(phi_intp))) .^ 2 - (sin(deg2rad(Phi_intp))) .^ 2) ./ (cos(deg2rad(phi_intp)));
u = imgaussfilt(u, 1);
T = Th_pad .* (Pres_intp / (pb / 100)) .^ kappa;
T0 = Th_pad .* (Pres0_intp / (pb / 100)) .^ kappa;
DTDH = zeros(size(Pres_intp));
DTDH0 = zeros(size(Pres_intp));
DTDH(:, 3:end - 2) = (T(:, 4:end - 1) - T(:, 2:end - 3)) ./ ...
    (R * T(:, 3:end - 2) / gr .* log(Pres_intp(:, 2:end - 3) ./ Pres_intp(:, 4:end - 1))) * 1000;
DTDH0(:, 3:end - 2) = (T0(:, 4:end - 1) - T0(:, 2:end - 3)) ./ ...
    (R * T0(:, 3:end - 2) / gr .* log(Pres0_intp(:, 2:end - 3) ./ Pres0_intp(:, 4:end - 1))) * 1000;

%% Plot
figure(1)
clf
set(gca, "Box", "on", "LineWidth", 1, "FontSize", 9)
hold on
contourf(phi_intp, Pres0_intp, Th_pad, ...
    300:0.5:311, "LineStyle", "none")
contour(phi_intp, Pres0_intp, Th_pad, ...
    300:1:311, "LineColor", "k", "LineWidth", 0.5)
hold off
cbar = colorbar("Box", "on", "LineWidth", 1, "FontSize", 9);
% colormap()
clim([300, 311])
ylabel(cbar, 'Potential temperature [K]', "FontSize", 11)
xlim([phis, phin])
ylim([600, 1000])
xticks(-30:10:30)
xticklabels({'30°S', '20°S', '10°S', '0°', '10°N', '20°N', '30°N'})
yticks([600, 750, 850, 900, 950, 1000])
xlabel('Latitude', "FontSize", 11)
ylabel('Pressure [hPa]', "FontSize", 11)
box on
set(gca, "Ydir", "reverse", "YScale", "log", "Layer", "top")
set(gcf, "Position", [10 10 440 240], "Color", "w")

if export
    export_fig(gcf, 'Figures/exp2_theta_initial', "-pdf", "-nocrop", "-nofontswap", "-painters")
end

figure(2)
clf
set(gca, "Box", "on", "LineWidth", 1, "FontSize", 9)
hold on
contourf(phi_intp, Pres_intp, Th_pad, ...
    300:0.5:311, "LineStyle", "none")
contour(phi_intp, Pres_intp, Th_pad, ...
    300:1:311, "LineColor", [0.5, 0.5, 0.5], "LineWidth", 0.5)
[C, h] = contour(phi_intp, Pres_intp, u, ...
    -1:0.1:1, "LineColor", "k", "LineWidth", 0.5);
clabel(C, h, -1:0.2:1)
hold off
cbar = colorbar("Box", "on", "LineWidth", 1, "FontSize", 9);
% colormap()
clim([300, 311])
ylabel(cbar, 'Potential temperature [K]', "FontSize", 11)
xlim([phis, phin])
ylim([600, 1000])
xticks(-30:10:30)
xticklabels({'30°S', '20°S', '10°S', '0°', '10°N', '20°N', '30°N'})
yticks([600, 750, 850, 900, 950, 1000])
xlabel('Latitude', "FontSize", 11)
ylabel('Pressure [hPa]', "FontSize", 11)
box on
set(gca, "Ydir", "reverse", "YScale", "log", "Layer", "top")
set(gcf, "Position", [10 250 440 240], "Color", "w")

if export
    export_fig(gcf, 'Figures/exp2_theta_final', "-pdf", "-nocrop", "-nofontswap", "-painters")
end

figure(3)
clf
set(gca, "Box", "on", "LineWidth", 1, "FontSize", 9)
hold on
plot(DTDH0(82, 3:end - 2), Pres0_intp(82, 3:end - 2), ...
    "Color", "k", "LineWidth", 0.5, "DisplayName", 'initial')
plot(DTDH(82, 3:end - 2), Pres_intp(82, 3:end - 2), ...
    "Color", "k", "LineWidth", 0.5, "Linestyle", "--", "DisplayName", 'final')
hold off
legend("show", "Box", "off", "FontSize", 11, "Location", "southeast")
xlim([-10, 5])
ylim([600, 1000])
yticks([600, 750, 850, 900, 950, 1000])
xlabel('Lapse rate [K/km]', "FontSize", 11)
ylabel('Pressure [hPa]', "FontSize", 11)
box on
set(gca, "Ydir", "reverse", "YScale", "log", "Layer", "top")
set(gcf, "Position", [450 10 160 240], "Color", "w")

if export
    export_fig(gcf, 'Figures/exp2_lapserate_0', "-pdf", "-nocrop", "-nofontswap", "-painters")
end

figure(4)
clf
set(gca, "Box", "on", "LineWidth", 1, "FontSize", 9)
hold on
plot(DTDH0(41, 3:end - 2), Pres0_intp(41, 3:end - 2), ...
    "Color", "k", "LineWidth", 0.5, "DisplayName", 'initial')
plot(DTDH(41, 3:end - 2), Pres_intp(41, 3:end - 2), ...
    "Color", "k", "LineWidth", 0.5, "Linestyle", "--", "DisplayName", 'final')
hold off
legend("show", "Box", "off", "FontSize", 11, "Location", "southeast")
xlim([-10, 5])
ylim([600, 1000])
yticks([600, 750, 850, 900, 950, 1000])
xlabel('Lapse rate [K/km]', "FontSize", 11)
ylabel('Pressure [hPa]', "FontSize", 11)
ylabel('Pressure [hPa]', "FontSize", 11)
box on
set(gca, "Ydir", "reverse", "YScale", "log", "Layer", "top")
set(gcf, "Position", [450 250 160 240], "Color", "w")

if export
    export_fig(gcf, 'Figures/exp2_lapserate_30', "-pdf", "-nocrop", "-nofontswap", "-painters")
end
