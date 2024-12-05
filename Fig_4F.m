%% function Fig_4F

clear
close all
clc

addpath 'helper_functions'

%% 1. Specify folder paths and load the structs
parent_dir = fullfile(pwd);
abs_path = fullfile(parent_dir, 'Visual_Weights/');

alpha_data = load(fullfile(abs_path, 'alpha_vector_smoothed.mat'));
alpha_data = alpha_data.alpha_vector_smoothed;

% alpha_data = load(fullfile(abs_path, 'alpha_vector_unsmoothed.mat'));
% alpha_data = alpha_data.alpha_vector_unsmoothed;

luminance_val = load(fullfile(abs_path, 'Luminance_cell.mat'));
luminance_val = luminance_val.Luminance_cell;
% exist(abs_path)


fishNames = {'Hope', 'Len', 'Doris', 'Finn', 'Ruby'}; % consistent with SICB
num_fish = length(fishNames);



colorMap = [51,160,44;
            201,108,255;
            255,80,11;
            32,81,178;
            179,0,0]/255;

all_lux = [];
all_alpha = [];



for i = 1 : num_fish

    alpha_data{i}
    all_lux = [all_lux, luminance_val{i}];
    all_alpha = [all_alpha; alpha_data{i}(:)];

end

%% Curve Fitting

X1 = log10(all_lux);
Y1 = all_alpha;
    
X1 = X1(:);
    
Y1 = Y1(:);

[xData, yData] = prepareCurveData( X1(:), Y1(:));

ft = fittype( 'a/(1+exp(-b*(x-c)))+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.33 2.275 0.6186 0.3275] + 0*0.1*[0.33 2.275 0.6186 0.3275].*randn(1,4);
    
 [fitresult1, gof1] = fit( xData, yData, ft, opts );



n_points = 100;
x_in = linspace(min(X1), max(X1), n_points);
y_in = feval(fitresult1, x_in);

ci = predint(fitresult1,x_in,0.95,'functional','on');
xx = [x_in,fliplr(x_in)];
yy = [ci(:,1);flipud(ci(:,2))]';


figure('Color','white')
set(gca,'LineWidth',1.5,'fontsize',14)
hold on

h(7) = fill(10.^(xx),yy,0.5*[1 1 1],'EdgeColor','none','FaceAlpha',0.5);

for i = 1 : num_fish

    
    h(i) = scatter((luminance_val{i}), (alpha_data{i}),80,colorMap(i, :),'filled','MarkerFaceAlpha',0.9);
end

h(6) = plot(10.^(x_in),y_in,'k--','LineWidth',2);

ipt = findchangepts(y_in,'MaxNumChanges',2,'Statistic','mean');

xlabel('Illumination, lx')
ylabel('Visual Weight, \alpha(\lambda)')
set(gca,'Xscale','log')
xlim([0.08 280])
ylim([0.2 0.7])
lux_ticks = [0.1, 1, 10, 100];
    
xticks(lux_ticks);
    
xticklabels(lux_ticks);
pos = get(gca,'Position');
height_cm = pos(4);
desired_length = 0.015; %cm
normalized_length = desired_length./height_cm;
set(gca,'TickLength', [normalized_length, 0.01],'TickDir','in')
legend(h,{'fish 1','fish 2','fish 3','fish 4','fish 5','sigmoidal fit','95% PI'},'Location','southeast','edgecolor','none','fontsize',16)

