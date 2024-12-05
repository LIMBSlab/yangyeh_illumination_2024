%% function Fig_2J

clear
close all
clc


close all
addpath 'helper_functions'



colorMap = [51,160,44
            201,108,255
            255,80,11
            32,81,178
            179,0,0]/255;


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
                res(i).luminances(il).rms(trial_idx) = rms_displacement.^2;
            end
        % end

        res(i).luminances(il).rms_mean = mean(res(i).luminances(il).rms);
        res(i).luminances(il).rms_std = std(res(i).luminances(il).rms, 'omitnan');

        % 05.18 NEW: standard error about mean
        res(i).luminances(il).rms_sem = std(res(i).luminances(il).rms, 'omitnan') / sqrt(num_trials);

    end
end


%% 3. Plotting starts here
num_body_pts = 12;
num_fish = 5;
target_pt = 12; % only look at tail
field_name = 'rmsMean';
% colorMap = cool(num_fish+1);

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

for i = 1:5

   

    fish_name = fishNames{i};

    % if i == 1 % Get rid of Hope lux 1, 3, 9
    %     lux = [0.4, 2, 3.5, 5.5, 7, 9.5, 15, 30, 60, 150, 210];
    % 
    % else
        lux = [res(i).luminances.lux];
    % end
    lux_valid{i} = lux;

    % Centered, then smoothed for x-variance values
    data = [res(i).luminances.rms_mean];
    % std_dev = [res(i).luminances.rms_std];
    std_sem = [res(i).luminances.rms_sem];

    data_smoothed_centered{i} = movmean(data - 1*(avg_rms(i) - mean_value_all), 3)';

    % Collect all data for Sigmoid fitting
    if i~=3
        all_lux = [all_lux, lux];
        all_data_pts_processed = [all_data_pts_processed; data_smoothed_centered{i}];
    end

end

    %% 4. Fit a sigmoid model
    X1 = log10(all_lux);
    Y1 = all_data_pts_processed;
    X1 = X1(:);
    Y1 = Y1(:);

 [fitresult1, gof1] = createSigmoidFit(field_name, X1, Y1);

    num_points = 500;
    x_in = linspace(min(X1), max(X1), num_points);

    % Evaluate the fitted model at the sample points
    y_in = feval(fitresult1, x_in);

    pi = predint(fitresult1,x_in,0.95,'functional','on');
    % pi = predint(fitted_model,x_in,0.95,'observation','on');
    xx = [x_in,fliplr(x_in)];
    yy = [pi(:,1);flipud(pi(:,2))]';



    figure('Color','white')
    set(gca,'LineWidth',1.5,'fontsize',14)
    hold on

    fill(10.^(xx),yy,0.5*[1 1 1],'EdgeColor','none','FaceAlpha',0.5);
   
    
    for i = 1 : num_fish
        h (i) = scatter((lux_valid{i}), (data_smoothed_centered{i}),80,colorMap(i, :),'filled','MarkerFaceAlpha',0.9);
    plot(lux_valid{i},data_smoothed_centered{i},'color',colorMap(i, :),'linewidth',2)
    end
     plot(10.^(x_in), y_in, 'Color', 'k', 'LineWidth', 2);
    xlabel('Illumination, lx')
    ylabel('Transverse Tail Movement Variance (cm^2)')
    set(gca,'Xscale','log')
    xlim([0.08 280])
    ylim([-0.2 2.2])
    lux_ticks = [0.1, 1, 10, 100];
    xticks(lux_ticks);
    xticklabels(lux_ticks);
    axis square
    legend(h,{'fish 1','fish 2','fish 3','fish 4','fish 5'},'FontSize',14)


return
    ipt = findchangepts(y_in,'MaxNumChanges',2,'Statistic','mean');
[10^x_in(ipt(1)),10^x_in(ipt(2))]

n_poly = 7;
for i = 1:n_poly
    p1{i} = polyfit(X1,Y1,i);
    % y{i} = polyval(p{i},x_in);
    SSE1(i) = sum((Y1-polyval(p1{i},X1)).^2);
    AIC_1(i) = length(X1)*log(SSE1(i)/length(X1)) + 2*(i+1+1);
    AICc_1(i) = AIC_1(i) + 2*(i+1)*(i+1+1)/(length(X1)-(i+1)-1);
   
    [~,r2adj_smooth(i)] = rsquared(Y1,polyval(p1{i},X1),i+1);

end

[r2_sigmoid_smooth,r2adj_sigmoid_smooth] = rsquared(Y1,feval(fitresult1, X1),4);
AIC_sigmoid_smooth = length(X1)*log(gof1.sse/length(X1)) + 2*(4+1);
AIC_sigmoid_smooth_c = AIC_sigmoid_smooth + 2*(4)*(4+1)/(length(X1)-(4)-1);


figure('Color','white')
set(gca,'LineWidth',1.5,'fontsize',14)
hold on
h(1) = plot(1:1:n_poly,r2adj_smooth,'b.-','MarkerSize',25,'LineWidth',2);

h(3) = plot(1:1:n_poly,r2adj_sigmoid_smooth*ones(1,n_poly),'k-','LineWidth',2);

xlabel('Degree of Polynomial')
ylabel('Adjusted R^2')
xlim([0.5 (n_poly+0.5)])
% ylim([0.55 0.85])


figure('Color','white')
set(gca,'LineWidth',1.5,'fontsize',14)
hold on
h(1) = plot(1:1:n_poly,AICc_1,'b.-','MarkerSize',25,'LineWidth',2);

h(3) = plot(1:1:n_poly,AIC_sigmoid_smooth_c*ones(1,n_poly),'k-','LineWidth',2);

xlabel('Degree of Polynomial')
ylabel('Small Sample Size Corrected AIC')
xlim([0.5 (n_poly+0.5)])
% ylim([0.55 0.85])

[min_val_R2,min_idx_R2] = min(r2adj_sigmoid_smooth - r2adj_smooth);
[r2adj_sigmoid_smooth,min_val_R2,min_idx_R2]
disp('--------------------------------------------------------------')
[min_val_AICc,min_idx_AICc] = min(AICc_1-AIC_sigmoid_smooth_c);
[AIC_sigmoid_smooth_c,min_val_AICc,min_idx_AICc]

return
%%
close all
clc
for i = 1:5

    
   
    clear X1 Y1 fitresult1 gof1
    X1 = log10(lux_valid{i});
    Y1 = data_smoothed_centered{i};
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

    if gof1.adjrsquare <= 0.89
        n_poly = length(Y1)-1;
        for k = 1:n_poly
            p1{k} = polyfit(X1,Y1,k);
            [~,r2adj_smooth_fish3(k)] = rsquared(Y1,polyval(p1{k},X1),k+1);

        end
    end
    figure('Color','white')
    set(gca,'LineWidth',1.5,'fontsize',14)
    hold on

    fill(10.^(xx),yy,0.5*[1 1 1],'EdgeColor','none','FaceAlpha',0.5);
   
    [i, gof1.adjrsquare]
     plot(lux_valid{i},data_smoothed_centered{i},'.-','color',colorMap(i, :),'linewidth',2,'MarkerSize',35)
    
     plot(10.^(x_in), y_in, 'Color', 'k', 'LineWidth', 2);
    
    xlabel('Illumination, lx')
    ylabel('Tail RMS (cm)')
    set(gca,'Xscale','log')
    xlim([0.08 280])
    ylim([0.2 1.6])
    lux_ticks = [0.1, 1, 10, 100];
    xticks(lux_ticks);
    xticklabels(lux_ticks);
    axis square
    title(fishNames{i},'FontWeight','normal')   

end




%%
 function [fitresult, gof] = createSigmoidFit(field_name, x, y)
%  Auto-generated by MATLAB on 12-Mar-2024 15:48:24

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
