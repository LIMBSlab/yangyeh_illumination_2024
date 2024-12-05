%% function Fig_S1A

clear
close all
clc

addpath 'helper_functions'

%% 1. Specify folder paths and load the structs
parent_dir = fullfile(pwd);
abs_path = fullfile(parent_dir, 'data_structures');


shuttle = load(fullfile(parent_dir, 'shuttle.mat'), 'shuttle').shuttle;
fishNames = {'Hope', 'Len', 'Doris', 'Finn', 'Ruby'}; % consistent with SICB
num_fish = 5;

Fs = 25;
TIME = 0:1/Fs:20-1/Fs;



all_fish = load(fullfile(abs_path, 'data_clean_head.mat'), 'all_fish').all_fish;
% Struct field names
x_field = 'fishX';
x_field_mean = 'fishXMean';



fish_name = fishNames{1};


% Set up colors
opacity = 0.25;
colors = [copper(7), (ones(7, 1)* opacity)];
color_idx = [1, 4, 6];
colors_mean = {'#2256A5', '#2C68F5', '#1FA9FF'}; 
colors_il = {'#543F36', '#805749', '#AC6C23'};

%% 2. Loop through and rotate all fish w.r.t. main body axis
fish_idx = queryStruct(all_fish, 'name', fish_name);
corr_val = nan(30,numel(all_fish(fish_idx).data));
il_range = 1:1:numel(all_fish(fish_idx).data); % 0.1, 2, 150 lux.
count = 1;

[Ps,freq] = single_side_specta(shuttle*100,1,Fs,0,[0 5]);
for il = il_range
% for il = [1, 14]
    % Each il gets a figure
   

    color_this_il = colors_mean{1};
    color_mean = colors_mean{1};
    

    % Plot the individual lines all in the same plot
    for tr_idx = 1 : numel(all_fish(fish_idx).data(il).(x_field))

        data = cell2mat(all_fish(fish_idx).data(il).(x_field)(tr_idx)); 
        data2 = fixSmallTL(data, 10)*100;
        
        [P1{il}(:,tr_idx),freq] = single_side_specta(data2,1,Fs,0,[0 5]);

       N = length(data2);  % Number of shuttle points
      fft_data2 = fft(data2);  % Compute the FFT
      frequencies = (0:N-1)*(Fs/N);  % Frequency axis in Hz
      fft_magnitude = abs(fft_data2/N);  % Magnitude of FFT (normalized)
      P2{il}(:,tr_idx) = fft_magnitude;


        corr_val(tr_idx,il) = corr(fixSmallTL(data, 10)*100,shuttle*100,'Type','Spearman');

    end
    P1_mean(:,il) = mean(P1{il},2,'omitmissing');
    P1_sem(:,il) = std(P1{il},[],2,'omitmissing')/sqrt(numel(all_fish(fish_idx).data(il).(x_field))-1);
    P2_mean(:,il) = mean(P2{il},2,'omitmissing');
    P2_sem(:,il) = std(P2{il},[],2,'omitmissing')/sqrt(numel(all_fish(fish_idx).data(il).(x_field))-1);

  
    
    ylim([-8 8]);


    color_shuttle = [1, 0, 0, 0.9]; % red
   
    data_mean = cell2mat(all_fish(fish_idx).data(il).(x_field_mean));

    trace_shuttle = plot(TIME, shuttle*100, 'color', color_shuttle, 'LineWidth', 3, 'DisplayName', 'Refuge');
    trace_data_mean = plot(TIME, fixSmallTL(data_mean, 10)*100, 'color', color_mean, 'LineWidth', 3, 'DisplayName', 'Average Data');



    corr_val_mean(il) = corr(shuttle*100,fixSmallTL(data_mean, 10)*100,'Type','Spearman');
    % legend([trace_shuttle, trace_data_mean, p1]);

   
    % Get title
    lux_measured(il) = all_fish(fish_idx).data(il).luxMeasured;
    num_trials(il) = numel(all_fish(fish_idx).data(il).fishX);
    
   

end
%%

%%



figure('Color','white');
    set(gca,'LineWidth',1.5,'FontSize',14)
    hold on
data_corr = smooth(corr_val_mean,5);
p = linspace(0,1,29);
colorMap = zeros(length(p),3);
colorMap(:,1) = p;



set(gca,'Xscale','log')
xlim([0.08 250])
ylim([0.75 0.95])
ylabel('Spearman Correlation Coefficient')
xlabel('Illumination (lx)')
pos = get(gca,'Position');
height_cm = pos(4);
desired_length = 0.015; %cm
normalized_length = desired_length./height_cm;
set(gca,'TickLength', [normalized_length, 0.01],'TickDir','in')
lux_ticks = [0.1, 1, 10, 100];
    xticks(lux_ticks);
    xticklabels(lux_ticks);


 X1 = log10(lux_measured(2:end));
 Y1 = data_corr(2:end);
[fitresult1, gof1] = createSigmoidFit(X1, Y1);
num_points = 500;
    x_in = linspace(log10(lux_measured(1)), log10(lux_measured(end)), num_points);

    % Evaluate the fitted model at the sample points
    y_in = feval(fitresult1, x_in);

    pi = predint(fitresult1,x_in,0.95,'functional','on');
    % pi = predint(fitted_model,x_in,0.95,'observation','on');
    xx = [x_in,fliplr(x_in)];
    yy = [pi(:,1);flipud(pi(:,2))]';

    h(2) = fill(10.^(xx),yy,0.5*[1 1 1],'EdgeColor','none','FaceAlpha',0.5);
    h(1) = plot(10.^(x_in), y_in, 'Color', 'k', 'LineWidth', 2);
    for i = 1:length(data_corr)
    scatter(lux_measured(i), data_corr(i),60,0*colorMap(2*i, :),'filled','MarkerFaceAlpha',0.9);
    end


ipt = findchangepts(y_in,'MaxNumChanges',2,'Statistic','mean');

[10^x_in(ipt(1)),10^x_in(ipt(2))]



%% Helper: Find struct by field name
function i = queryStruct(struct, fieldName, query)
for i = 1:numel(struct)
    if isfield(struct(i), fieldName) && isequal(struct(i).(fieldName), query)
        return;
    end
end
end

%% Helper: fixSmallTL
% Updated 09.19.2022 by Joy Yeh
% Fill up small tracking losses after applying a 5Hz cut-off butterworth
% filter
%
% Params:
% data: the time-domain x value passed in. Might have small tracking loss
% window: the movmedian window size (10 is good)
%
% Returns:
% fixedData: filled gap with movmedian and applied a butterworth filter.
function [fixedData] = fixSmallTL(data, windowSize)
filled = fillmissing(data, 'movmedian', windowSize);

if sum(isnan(filled)) > 0
    fixedData = filled;
    return;
else
    % butterworth
    fs = 25;
    fc = 5; % 5 Hz cutoff
    Wn = fc /(fs / 2); % Cut-off for discrete-time filter
    [b,a] = butter(2, Wn);
    fixedData = filtfilt(b,a, filled);
end
end


function [fitresult, gof] = createSigmoidFit(x, y)

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a/(1+exp(-b*(x-c)))+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';


    opts.StartPoint = [0.0905 10.29 1.002 0.8011];


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end


