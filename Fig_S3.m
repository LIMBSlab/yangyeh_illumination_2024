%% function Fig_S3

clear
close all
clc


close all
addpath 'helper_functions'

%% 1. Specify folder paths
parent_dir = fullfile(pwd);
abs_path = fullfile(parent_dir, 'data_structures');


%% 2. Initial setup
fishNames = {'Fish 1', 'Fish 2', 'Fish 3', 'Fish 4', 'Fish 5'}; % consistent with SICB
numIls = [14, 9, 11, 9, 9];
numFish = 5;

shuttle = load('shuttle.mat');
shuttle = shuttle.shuttle;

c = copper(10);
colorMap = [c(4, :); c(6, :); c(7, :); c(8, :); c(9, :)];
gray = [0.7, 0.7, 0.7];

% Struct field names
gainField = 'gmGain';
phaseField = 'gmPhase';

all_fish = load(fullfile(abs_path, 'data_clean_head.mat'), 'all_fish').all_fish;

yLimG = [0 10];
yLimP = [-200 -22.5];

% Locate the frequency peaks
k = [2, 3, 5, 7, 11, 13, 19, 23, 29, 31, 37, 41];
freq_data = k * 0.05;

%% 2. Plotting starts here
colorMap = {'#000000', '#112d80', '#234099', '#2d50b4', '#5070c7', ...
    '#91a6e2', '#acc0fa', '#cbd0ee', '#ececec'};


for fish_idx = 1:5
    fish_name = fishNames{fish_idx};

    figure('Color','white');
    subplot(211)
    set(gca,'LineWidth',1,'fontsize',12)
    hold on
    % f.Position = [100 100 300 550];
    il = 0;

    num_il_levels = numel(all_fish(fish_idx).data);
    colorMap = magma(num_il_levels+1);

    plotName = [fish_name, ' Closed-Loop Frequency Responses'];

    % % Gain
    % lineWidth = 1.9;
    % h1 = axes('position',[0.2 0.56 0.76 0.4]);
    % hold on

    for il = 1 : num_il_levels
        data = all_fish(fish_idx).data(il).(gainField);
        c = colorMap(il, :);
        semilogx(freq_data, smooth(data), 'color', c, 'LineWidth',2);
    end

 
    xlim([0 2.1]);

    % set(h1, 'XTick', [0.1, 1]);
    % h1.XAxis.FontSize = labelFontSize;

    set(gca,'xScale','log');
    set(gca,'yScale','log');


    ylim([0.08, 1.5]);


    ylabel('Bode Gain (cm/cm)', 'FontSize', 12)
    title(plotName);

    
    subplot(212)
    set(gca,'LineWidth',1,'fontsize',12)
    hold on
    
    for il = 1 : num_il_levels
        data = all_fish(fish_idx).data(il).(phaseField);
        c = colorMap(il, :);
        semilogx(freq_data, smooth(data), 'color', c, 'LineWidth', 2);
    end


    xlim([0 2.1]);
    set(gca,'xScale','log');


    xlabel('Frequency (Hz)', 'FontSize', 12);
  
    ylim([-210, 0]);
    yticks([-210 -180 -150 -120 -90 -60 -30 0 30]);

    
    ylabel('Bode Phase (deg)', 'FontSize', 12);
    % colorbar
    colormap(colorMap(1:end-1,:))
    
end