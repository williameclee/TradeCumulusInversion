figure(2)
clf
subplot(2, 2, 1)
contourf(Phiu(2:2:end-1), theta_pad(2:end - 1), Pres(2:2:end-1, 2:end - 1).' * pb / 100, 100:50:1000)
xlim([phis, phin])
ylim([thetab, thetat])
colorbar

subplot(2, 2, 2)
contourf(Phiu(2:2:end - 1), theta_pad(2:end - 1), X(2:2:end - 1, 2:end - 1).' * cs / 1000)
xlim([phis, phin])
ylim([thetab, thetat])
colorbar

subplot(2, 2, 3)
contourf(Phiu(1:2:end), theta_pad(2:end - 1), X(1:2:end, 2:end - 1).')
xlim([phis, phin])
ylim([thetab, thetat])
colorbar

% subplot(2, 2, 4)
% contourf(Phiu(3:2:end-2), theta_pad(2:end - 1), Sigma(3:2:end-2, 2:end - 1).' * sigma0)
% xlim([phis, phin])
% ylim([thetab, thetat])
% colorbar