%% function Fig3C

clear
close all
clc

addpath 'helper_functions'

%% 1. Specify folder paths and load the structs
parent_dir = fullfile(pwd);
abs_path = fullfile(parent_dir, 'data_structures');

% exist(abs_path)

all_fish = load(fullfile(abs_path, 'data_clean_head.mat'), 'all_fish').all_fish;
fishNames = {'Hope', 'Len', 'Doris', 'Finn', 'Ruby'}; % consistent with SICB
fishNames_str = {'Fish 1', 'Fish 2', 'Fish 3', 'Fish 4', 'Fish 5'}; % consistent with SICB
p2m = 0.0004;
  
num_body_pts = 12;
num_fish = 5;
target_pt = 12; % only look at tail
field_name = 'gmPhase';


colorMap = [51,160,44
            201,108,255
            255,80,11
            32,81,178
            179,0,0]/255;

%% 2. Extract data
all_lux = [];
all_data_pts = [];
all_data_pts_processed_1 = [];
all_data_pts_processed_2 = [];
avg_phase = zeros(5, 1);

data_cell = cell(5, 14);
for i = 1 :num_fish
    num_ils = numel(all_fish(i).data);

    for il = 1 : num_ils
        phases = [all_fish(i).data(il).gmPhase];
        phases_mean = mean(phases(10:end)); % Frequencies 1.55 - 2.05 Hz
        data_cell{i, il} = phases_mean;
    end   

    data_smoothed = smooth([data_cell{i, :}],3)';
    avg_phase(i) = mean(data_smoothed, 'omitnan');
    all_data_pts = [all_data_pts, data_smoothed];
end
all_data_mean = mean(all_data_pts);


lux_level = [0.1	0.4	1	2	3.5	5.5	7	9.5	12	15	30	60	150	210];






for i = 1 : num_fish
    fish_name = fishNames{i};

    % Calculate centered, smoothed, and population mean-subtracted phase
    lux{i} = [all_fish(i).data.luxMeasured];

    data_1 = smooth([data_cell{i, :}]);
    data_centered_1{i} = (data_1 - 1*(avg_phase(i) - 1*all_data_mean));

    data_2 = smooth([data_cell{i, :}],1);
    data_centered_2{i} = (data_2 - 1*(avg_phase(i) - 1*all_data_mean));

    % Collect all data for Sigmoid fitting
    all_lux = [all_lux, lux{i}];
    all_data_pts_processed_1 = [all_data_pts_processed_1; data_centered_1{i}];
    all_data_pts_processed_2 = [all_data_pts_processed_2; data_centered_2{i}];

    
end






%% Curve Fitting
X1 = log10(all_lux);
Y1 = all_data_pts_processed_1;
X1 = X1(:);
Y1 = Y1(:);

X2 = log10(all_lux);
Y2 = all_data_pts_processed_2;
X2 = X2(:);
Y2 = Y2(:);






[xData1, yData1] = prepareCurveData( X1(:), Y1(:));
[xData2, yData2] = prepareCurveData( X2(:), Y2(:));

ft = fittype( 'a/(1+exp(-b*(x-c)))+d', 'independent', 'x', 'dependent', 'y' );
opts1 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts1.Display = 'Off';

opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
% for i = 1:10
    opts1.StartPoint = [34 1.21 2 -160] + 0*0.1*[34 1.21 2 -160].*randn(1,4);
    opts2.StartPoint = [33.8818    1.1134    2.2037 -171] + 0*0.1*[34 1.21 2 -160].*randn(1,4);
    [fitresult1, gof1] = fit( xData1, yData1, ft, opts1 );
    % [fitresult2, gof2] = fit( xData2, yData2, ft, opts2 );
    load("fig4B_sigmoidfit.mat")

%         res = Y-fitresult(X);
% [h,pValue,stat,cValue] = archtest(res)



    
    n_points = 100;
    x_in1 = linspace(min(X1), max(X1), n_points);
    y_in1 = feval(fitresult1, x_in1);

    x_in2 = linspace(min(X2), max(X2), n_points);
    y_in2 = feval(fitresult2, x_in2);


    
% end
ci_1 = predint(fitresult1,x_in1,0.95,'functional','on');
xx1 = [x_in1,fliplr(x_in1)];
yy1 = [ci_1(:,1);flipud(ci_1(:,2))]';

ci_2 = predint(fitresult2,x_in2,0.95,'functional','on');
xx2 = [x_in2,fliplr(x_in2)];
yy2 = [ci_2(:,1);flipud(ci_2(:,2))]';



figure('Color','white')
set(gca,'LineWidth',1.5,'fontsize',14)
hold on

h(7) = fill(10.^(xx1),yy1,0.5*[1 1 1],'EdgeColor','none','FaceAlpha',0.5);
for i = 1 : num_fish
    h(i) = scatter(lux{i}, (data_centered_1{i}),80,colorMap(i, :),'o','filled','MarkerFaceAlpha',0.9,'LineWidth',1);
end


h(6) = plot(10.^(x_in1),y_in1,'k--','LineWidth',2);

xlabel('Illumination, lx')
ylabel('High Frequency Tracking Phase (deg)')
set(gca,'Xscale','log')
xlim([0.08 280])

% set(gca,'Position', 'centimeters')
pos = get(gca,'Position');
height_cm = pos(4);
desired_length = 0.015; %cm
normalized_length = desired_length./height_cm;
set(gca,'TickLength', [normalized_length, 0.01],'TickDir','in')
lux_ticks = [0.1, 1, 10, 100];
    xticks(lux_ticks);
    xticklabels(lux_ticks);
legend(h,{'fish 1','fish 2','fish 3','fish 4','fish 5','sigmoidal fit','95% PI'},'Location','southeast','edgecolor','none','fontsize',16)

