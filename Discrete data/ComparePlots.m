
close all
plot(Bilinear0_2, 'linewidth', 1.5)
hold on
plot(Bilinear0_52, 'linewidth', 1)
plot(Bilinear0_02, 'linewidth', 1.5)
plot(reference_signal, 'linewidth', 1.5)
legend('T_s = 0.2', 'T_s = 0.52', 'T_s = 0.02', 'Reference signal')
xlim([0 10])
title('Bilinear Discretisation Step Response')
ylabel('Amplitude')
% 
% plot(Euler0_2, 'linewidth', 1.5)
% hold on
% plot(Euler0_52, 'linewidth', 1)
% plot(ans, 'linewidth', 1.5)
% plot(reference_signal, 'linewidth', 1.5)
% legend('T_s = 0.2', 'T_s = 0.52', 'T_s = 0.02', 'Reference signal')
% xlim([0 10])
% title('Euler Discretisation Step Response')
% ylabel('Amplitude')

% plot(ZOH0_2, 'linewidth', 1.5)
% hold on
% plot(ZOH0_46, 'linewidth', 1)
% plot(ZOH0_02, 'linewidth', 1.5)
% plot(reference_signal, 'linewidth', 1.5)
% legend('T_s = 0.2', 'T_s = 0.46', 'T_s = 0.02', 'Reference signal')
% xlim([0 10])
% title('Zero Order Hold Discretisation Step Response')
% ylabel('Amplitude')
