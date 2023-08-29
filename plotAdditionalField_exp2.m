%% Plot
CustomColour
figure(1)
clf
set(gca, "Box", "on", "LineWidth", 1, "FontSize", 9)
hold on
contourf(phi_intp(:,3:end-2), Pres_intp(:,3:end-2), DTDH0(:,3:end-2), ...
    -10:0.2:4, "LineStyle", "none")
[C, h] = contour(phi_intp(:,3:end-2), Pres_intp(:,3:end-2), DTDH0(:,3:end-2), ...
    -8:4:4, "LineColor", "k", "LineWidth", 0.5);
clabel(C, h, -8:4:4)
hold off
cbar = colorbar("Box", "on", "LineWidth", 1, "FontSize", 9);
colormap(ccmap.seismic)
clim([-10, 10])
ylabel(cbar, 'Lapse rate [K/km]', "FontSize", 11)
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
export_fig(gcf, 'Figures/exp2_lapserate_initial', "-pdf", "-nocrop", "-nofontswap", "-painters")

figure(2)
clf
set(gca, "Box", "on", "LineWidth", 1, "FontSize", 9)
hold on
contourf(phi_intp(:,3:end-2), Pres_intp(:,3:end-2), DTDH(:,3:end-2), ...
    -10:0.2:4, "LineStyle", "none")
[C, h] = contour(phi_intp(:,3:end-2), Pres_intp(:,3:end-2), DTDH(:,3:end-2), ...
    -8:4:4, "LineColor", "k", "LineWidth", 0.5);
clabel(C, h, -8:4:4)
hold off
cbar = colorbar("Box", "on", "LineWidth", 1, "FontSize", 9);
colormap(ccmap.seismic)
clim([-10, 10])
ylabel(cbar, 'Lapse rate [K/km]', "FontSize", 11)
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
export_fig(gcf, 'Figures/exp2_lapserate_final', "-pdf", "-nocrop", "-nofontswap", "-painters")

plot_Phi = [0,10,20,30,40];
plot_Phiid = [82,96,110,122,134];
fig_Num = 3;

for fig_num = 1:5
    figure(fig_num+2)
    clf
    set(gca, "Box", "on", "LineWidth", 1, "FontSize", 9)
    hold on
    plot(DTDH0(plot_Phiid(fig_num), 3:end - 2), Pres0_intp(plot_Phiid(fig_num), 3:end - 2), ...
        "Color", cc(4).blue, "LineWidth", 0.5, "DisplayName", 'initial')
    plot(DTDH(plot_Phiid(fig_num), 3:end - 2), Pres_intp(plot_Phiid(fig_num), 3:end - 2), ...
        "Color", cc(4).red, "LineWidth", 0.5, "DisplayName", 'final')
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
    export_fig(gcf, ['Figures/exp2_lapserate_',num2str(plot_Phi(fig_num),2)], "-pdf", "-nocrop", "-nofontswap", "-painters")
end