load SE_CP.mat
load SV_CP.mat

% Frequencies
f = 0.05 * [2 3 5 7 11 13 19 23 29 31 37 41];

% Plot the S_V(jw)C(jw)P(jw) and S_E(jw)C(jw)P(jw)
figure,
subplot 211
hold on
semilogx(f,smooth(abs(SE_CP)),'color','#369DAB','LineWidth',2);
semilogx(f,smooth(abs(SV_CP)),'color','#FFC000','LineWidth',2);
set(gca,'xScale','log');
set(gca,'yScale','log');
legend('G_E(s)','G_V(s)')
ylabel('Bode Gain (cm/cm)')
title('Fig 4H')


subplot 212
hold on
semilogx(f,smooth(rad2deg(unwrap(angle(SE_CP)))),'color','#369DAB','LineWidth',2);
semilogx(f,smooth(rad2deg(unwrap(angle(SV_CP)))),'color','#FFC000','LineWidth',2);
set(gca,'xScale','log');
xlabel('Frequency (Hz)')
ylabel('Bode Phase (deg)')