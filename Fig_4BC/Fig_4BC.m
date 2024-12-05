% load data from fish 1:
load data_clean_head_with_GM.mat

% load fitted parameters to forward-loop Cascade SCP
load Best_parameter_CP_Hope_Smooth_2.mat

struct = all_fish(1).data;
parameter = Best_parameter_CP_Hope_Smooth_2;

f = 0.05 * [2 3 5 7 11 13 19 23 29 31 37 41];
w = 2 * pi * f;


CP_data = cell(1,2);
flag = 1;
for idx = [1,14]
CP_data{flag} = struct(idx).CPMean;
flag = flag + 1;
end

% Lowest illumination (0.1 lx) model
CP_model_predict{1} = (parameter(14,1)*(1i*w) + parameter(14,2))./((1i*w).*((1i*w) + parameter(14,3))) .* exp(-parameter(1,4)*(1i*w));

% Highest illumination (210 lx) model
CP_model_predict{2} = (parameter(14,1)*(1i*w) + parameter(14,2))./((1i*w).*((1i*w) + parameter(14,3))) .* exp(-parameter(14,4)*(1i*w));
% Note that here only the delay is different, the rest parameters are the
% same.

fig = figure();
fig.Position = [100 100 300 550];

subplot 211
hold on
grid on

semilogx(f,smooth(abs(CP_data{2})),'square','color','#fd9567','LineWidth',2,'MarkerFaceColor','#fd9567');

semilogx(f,abs(CP_model_predict{2}),'color','#fd9567','LineWidth',2);


set(gca,'xScale','log');
set(gca,'yScale','log');
axis([0.1 2.1 0.2 3]);
legend('210 lx (Exp)','210 lx (Model)')
ylabel('Bode Gain (cm/cm)')
title('Fig 4B')


subplot 212
hold on
grid on
    semilogx(f,smooth(rad2deg(unwrap(angle(CP_data{2})))),'square','color','#fd9567','LineWidth',2, 'MarkerFaceColor','#fd9567');
    semilogx(f,(rad2deg(unwrap(angle(CP_model_predict{2})))),'color','#fd9567','LineWidth',2);
set(gca,'xScale','log');
axis([0.1 2.1 -200 0]);
xlabel('Frequency (Hz)')
ylabel('Bode Phase (deg)')



%%
fig = figure();
fig.Position = [100 100 300 550];

subplot 211
hold on
grid on

semilogx(f,smooth(abs(CP_data{1})),'ko','LineWidth',2,'MarkerFaceColor','#000000');

semilogx(f,abs(CP_model_predict{1}),'k','LineWidth',2);


set(gca,'xScale','log');
set(gca,'yScale','log');
axis([0.1 2.1 0.2 3]);
legend('0.1 lx (Exp)','0.1 lx (Model)')
ylabel('Bode Gain (cm/cm)')
title('Fig 4C')


subplot 212
hold on
grid on
    semilogx(f,smooth(rad2deg(unwrap(angle(CP_data{1})))),'ko','LineWidth',2,'MarkerFaceColor','#000000');
    semilogx(f,(rad2deg(unwrap(angle(CP_model_predict{1})))),'k','LineWidth',2);
set(gca,'xScale','log');
axis([0.1 2.1 -200 0]);
xlabel('Frequency (Hz)')
ylabel('Bode Phase (deg)')

