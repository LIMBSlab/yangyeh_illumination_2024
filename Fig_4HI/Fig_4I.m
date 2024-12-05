load SE_CP.mat
load SV_CP.mat

load data_clean_head_with_GM.mat

 

% Frequencies
f = 0.05 * [2 3 5 7 11 13 19 23 29 31 37 41];

%% Plot the frequency domain tracking error
figure,
hold on
semilogx(f,smooth(abs(1./(1 + all_fish(1).data(1).CPMean))),'color','#333333','LineWidth',2);
semilogx(f,smooth(abs(1./(1 + all_fish(1).data(14).CPMean))),'color','#FD9567','LineWidth',2);
semilogx(f,smooth(abs(1./(1 + SE_CP))),'color','#369DAB','LineWidth',2);
semilogx(f,smooth(abs(1./(1 + SV_CP))),'color','#FFC000','LineWidth',2);
yline(1,'--');
legend('0.1 lx (Exp)','210 lx (Exp)','From G_E(s)','From G_V(s)','Location','southeast');
set(gca,'xScale','log');
xlabel('Frequency (Hz)');
ylabel('Tracking Error');
title('Fig 4I')