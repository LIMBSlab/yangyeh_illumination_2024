%% function Fig_S1D


clear
close all
clc

format long  g
close all
addpath 'helper_functions'

%% 1. Specify folder paths
parent_dir = fullfile(pwd);
abs_path = fullfile(parent_dir, 'data_structures');

%% 2. Load the full body clean struct
all_fish = load(fullfile(abs_path, 'data_clean_body.mat'), 'all_fish').all_fish;
raw = load(fullfile(abs_path, 'result_tail_positions.mat'), 'raw').raw;
res = raw;

fishNames = {'Hope', 'Len', 'Doris', 'Finn', 'Ruby'}; % consistent with SICB
p2m = 0.0004;
num_frames = 500;

%% 2. Populate data calculations for RMS

for i = 1 : 5
    num_ils = numel(all_fish(i).luminance);
    for il = 1: num_ils

        res(i).luminances(il).rms = [];

        num_trials = numel(res(i).luminances(il).x_tail);
        num_trials_comb{i}(il) = num_trials;
        % if num_trials < 4
        %     continue;
        % else
            for trial_idx = 1 : num_trials
                y = cell2mat(res(i).luminances(il).y_tail(trial_idx));
                rms_displacement = rms((y - mean(y))*p2m*100,'omitnan');
                res(i).luminances(il).rms(trial_idx) = rms_displacement;
            end
        % end

        res(i).luminances(il).rms_mean = mean(res(i).luminances(il).rms.^2);
        res(i).luminances(il).rms_std = std(res(i).luminances(il).rms.^2, 'omitnan');

        % 05.18 NEW: standard error about mean
        res(i).luminances(il).rms_sem = std(res(i).luminances(il).rms.^2, 'omitnan') / sqrt(num_trials);

    end
end


%% 3. Plotting starts here
num_body_pts = 12;
num_fish = 5;
target_pt = 12; % only look at tail
field_name = 'rmsMean';
colorMap = cool(num_fish+1);

%% 2. Gather data
all_lux = [];
all_data_pts = [];
all_data_pts_processed = [];
avg_rms = zeros(5, 1);

for i = 1:num_fish
    num_ils = numel(res(i).luminances);

    rms_values = [res(i).luminances.rms_mean];

    % Calculate the average of the 'rms' field values
    avg_rms(i) = mean(rms_values, 'omitnan');

    all_data_pts = [all_data_pts, rms_values];

end

mean_value_all = mean(all_data_pts);

for i = [1]

    


    fish_name = fishNames{i};

    % if i == 1 % Get rid of Hope lux 1, 3, 9
    %     lux = [0.4, 2, 3.5, 5.5, 7, 9.5, 15, 30, 60, 150, 210];
    % else
        lux = [res(i).luminances.lux]
    % end

    % Centered, then smoothed for x-variance values
    data = [res(i).luminances.rms_mean];
    % std_dev = [res(i).luminances.rms_std];
    std_sem = [res(i).luminances.rms_sem];

    data_smoothed_centered = movmean(data - 0*(avg_rms(i) - mean_value_all), 3)';

    % Collect all data for Sigmoid fitting
    all_lux = [all_lux, lux];
    all_data_pts_processed = [all_data_pts_processed; data_smoothed_centered];

    

    lux_valid = lux;
    lux_valid([1 3 9]) = [];
    data_smoothed_centered_valid = data_smoothed_centered;
    data_smoothed_centered_valid([1 3 9]) = [];
    std_sem_valid = std_sem;
    std_sem_valid([1 3 9]) = [];
    xx2 = [lux_valid,fliplr(lux_valid)];
    y1 = data_smoothed_centered_valid + smooth(std_sem_valid,3);
    y2 = data_smoothed_centered_valid - smooth(std_sem_valid,3);
    yy2 = [y1(:);flipud(y2(:))]';
    

end

 X1 = log10(lux_valid);
 Y1 = data_smoothed_centered_valid;
 X1 = X1(:);
 Y1 = Y1(:);

  [fitresult1, gof1] = createSigmoidFit(field_name, X1, Y1);

num_points = 500;
    x_in = linspace(min(X1), max(X1), num_points);

    % Evaluate the fitted model at the sample points
    y_in = feval(fitresult1, x_in);

    pi = predint(fitresult1,x_in,0.95,'functional','on');
    xx = [x_in,fliplr(x_in)];
    yy = [pi(:,1);flipud(pi(:,2))]';
%%

    

    %%
x_in = linspace(log10(lux(1)), log10(lux(end)), num_points);

    % Evaluate the fitted model at the sample points
    y_in = feval(fitresult1, x_in);

    pi = predint(fitresult1,x_in,0.95,'functional','on');
    xx = [x_in,fliplr(x_in)];
    yy = [pi(:,1);flipud(pi(:,2))]';

    figure('color','white')
    set(gca,'linewidth',1,'fontsize',12)
    hold on
    fill(10.^(xx),yy,0.5*[1 1 1],'EdgeColor','none','FaceAlpha',0.5);
    plot(lux, data_smoothed_centered,'.','MarkerSize',25)
    plot(10.^(x_in), y_in, 'Color', 'k', 'LineWidth', 2);
    xlabel('Illumination, lx')
    ylabel('Transverse Tail Movement Variance (cm^2)')
    set(gca,'Xscale','log')
    xlim([0.08 280])
    % ylim([0 1.6])
    lux_ticks = [0.1, 1, 10, 100];
    xticks(lux_ticks);
    xticklabels(lux_ticks);


%%


function [fitresult, gof] = createSigmoidFit(field_name, x, y)


[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a/(1+exp(-b*(x-c)))+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

% if (field_name == "varX_mean")
    opts.StartPoint = [0.473288848902729 0.351659507062997 0.830828627896291 0.585264091152724];
% elseif (field_name == "varY_mean")
%     opts.StartPoint = [0.0357116785741896 0.849129305868777 0.933993247757551 0.678735154857773];
% end

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end