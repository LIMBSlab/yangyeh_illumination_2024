% Script: Time_Domain_Clean_Up_Mahal_RMS.m
% Add validity tags to fish structs
% Use Complex Mahalanobis Z-scores and TD x-pos RMS Error Z-scores.

% Updated 11/15/2022

% 1. Load fish struct and calculate FFT
clear all;
close all;

fishNameList = {'hope', 'len', 'doris', 'finn', 'ruby'};
N_FISH = 5;

% Prepare for plotting

parent_dir = fullfile(pwd);
   
abs_path = fullfile(parent_dir, 'data_structures\');
    load([abs_path, '\pre_cleaned_original_structs\shuttle.mat']);
    savePath = fullfile(parent_dir, 'figures\');
   

% IMPORTANT: Controls the selection of outlier trials
% If there are more than 2 peaks with mahal z-scores > 1.5, also give a
% warning
zThreshAccumulate = 1.5;
multiPeakThresh = 3;

zThreshold = 2.5; % A freq peak is considered "invalid" if the z-score here > 2.2
maxAllowedOutlierPeaks = 1; % If ANY peak has a z-score > zThreshold, take it out
zRmsThreshold = 2.2; % Z-score threshold for rms values

% Specity freq values
k = [2, 3, 5, 7, 11, 13, 19, 23, 29, 31, 37, 41];
FREQ_AMPLITUDE_SCALAR = 0.05;
freqData = k * FREQ_AMPLITUDE_SCALAR;

% All the calculations
for idx = 4 %[1, 3, 4, 5]
    fishName = fishNameList{idx};

    % [NEW] the following line replaced "BigStruct"
    % load(['../Fish_Structs/', fishName, 'BigStruct.mat']);
    % load([pwd, '\11-15_Mahal_RMS_Tuned\Original\', fishName, 'Original.mat']);
    % load('../shuttle.mat');
    
    load([abs_path, '\pre_cleaned_original_structs\', fishName, 'BigStruct.mat']);
    

    numConditions = size(fish, 2);

    for il = 1 %1 : numConditions

        thisCell = fish(il).fishXClean;
        numTrials = size(thisCell, 2);

        [fftValsAll, fftMagsAll] = calculateFFTbyTrials(thisCell);
        [mahalDistances, Zmahal] = calculateMahalZScores(fftValsAll);
        [rmsErrors, Zrms] = calculateRmsZScores(shuttle, thisCell);

        % [NEW] IF more than 3 peaks have mahal z-scores greater than 1.5,
        % it's also not ok

        multiPeakWarning = zeros(numTrials, 1);
        for p = 1 : numTrials
            multiPeakWarning(p) = sum(abs(Zmahal(p, :)) > zThreshAccumulate);
        end

        % 4. Toggle the validity (N x 12)
        accumulateOK = multiPeakWarning < multiPeakThresh;
        warning = abs(Zmahal) > zThreshold;
        rmsOK =  abs(Zrms) < zRmsThreshold;

        % [NEW] Added 'accumulateOK' as a metric for outlier detection
        % validity = and((sum(warning, 2) < maxAllowedOutlierPeaks + 1),
        % rmsOK);

        validity = and((sum(warning, 2) < maxAllowedOutlierPeaks + 1), ...
            and(rmsOK, accumulateOK));
        numOutliers = numTrials - sum(validity);

        % 5. Save validity matrix back into the fish struct
        fish(il).xPosValidity = validity;
        fish(il).numOutliers = numOutliers;

        % 6. Get the time-domain data for the second-round cleaning
        fish(il).xClean02 = {};
        fish(il).xClean02Tr = {};
        fish(il).xClean02Rep = {};

        for j = 1 : numTrials
            if fish(il).xPosValidity(j) == 1
                fish(il).xClean02 = [fish(il).xClean02, fish(il).fishXClean{j}];
                fish(il).xClean02Tr = [fish(il).xClean02Tr, fish(il).fishXCleanTr{j}];
                fish(il).xClean02Rep = [fish(il).xClean02Rep, fish(il).fishXCleanRep{j}];
            end
        end

        % 7. Add the new average t-d field.
        XMean = cell(1, 1);
        currXTrials = fish(il).xClean02;
        XMean{1} = nanmean(cell2mat(currXTrials), 2);
        fish(il).xClean02Mean = XMean;

        % 8. Calculate the new freq response values
        xMeanData = fish(il).xClean02Mean;
        % disp(['Tick', ' ', num2str(il)]);
        [GM, gain, phase, cp, cpGain, cpPhase] = calculateFreqResponseValues(shuttle, xMeanData);

        fish(il).GM02 = GM;
        fish(il).gainMean02 = gain;
        fish(il).phaseMean02 = phase;

        fish(il).cp02 = cp;
        fish(il).cpGain02 = cpGain;
        fish(il).cpPhase02 = cpPhase;

        % Forced fix: rename the cpClean02 to cp here
        % Bug fix on 11/30/2022
        fish(il).cpClean02 = cp;

        fish(il).fftMagsAll = fftMagsAll;

    end


  
    zString = strrep(num2str(zThreshold), '.', '_');
    zRmsString = strrep(num2str(zRmsThreshold), '.', '_');
    zAccString = strrep(num2str(zThreshAccumulate), '.', '_');
    
    % thisFigFolder = [savePath, 'z_', zString, '_lim_', num2str(maxAllowedOutlierPeaks) '_zRMS_', ...
    %     zRmsString, '_acc_thresh_', zAccString, '_enum_', num2str(multiPeakThresh), '\'];
    % mkdir(thisFigFolder);

    for il = 1 %1 : numConditions
        numTrials = size(fish(il).fishXClean, 2);
        plotName = [fishName, ' Il = ', num2str(il), ' Fish X-Positions. ', num2str(numTrials), ' trials total.'];
        yLim = [-0.1 0.1];
        TIME = 0:0.04:20-0.04;

        %% Plotting starts here

        % Plot FFT
        plotFFT = 0;
        if (plotFFT == 1)
         figure;

         % Set a larger figure window for better visualization
         set(gcf, 'Position', [100, 100, 800, 400]);  % [left, bottom, width, height]

         for trNo = 1:numTrials
             hold on
             plot(freqData, fish(il).fftMagsAll(trNo, :), 'Color', [0.2 0.2 0.6 0.4], 'LineWidth', 1.5);
         end

         outlierIdx = find(fish(il).xPosValidity == 0);
         numOutliers = length(outlierIdx);

         legendInfo = {};
         for j = 1:numOutliers
             % Highlight outliers in red
             plot(freqData, fish(il).fftMagsAll(outlierIdx(j), :), 'Color', 'r', 'LineWidth', 2);
             legendInfo{j} = ['Trial # ', num2str(outlierIdx(j))]; % Customize legend info for outliers
         end

         % Add x and y axis labels
         xlabel('Frequency (Hz)');  % Customize as needed
         ylabel('Magnitude (cm)');  % Customize as needed

         % Adjust the x-axis tick labels to avoid cluttering
         xticks(freqData);
         xticklabels(freqData);
         set(gca, 'XTickLabelRotation', 45);  % Rotate x-axis labels to 45 degrees for better readability

         % Set appropriate axis limits if needed
         xlim([0, 2.10]);

         % Optional: Add grid for better visualization
        
         % Save the figure as both PNG and PDF
         saveas(gcf, [savePath, fishName, '_il_', num2str(il), '_FFT.png']);
         saveas(gcf, [savePath, fishName, '_il_', num2str(il), '_FFT.pdf']);

         % Set title for the plot
         title(plotName);
     end

        plotTimeDomain= 1;
        if plotTimeDomain == 1
            fig = figure();
            set(fig,'defaultLegendAutoUpdate','off');
            hold on

            % Highlight the outliers
            outlierIdx = find(fish(il).xPosValidity == 0);
            numOutliers = length(outlierIdx);

            legendInfo = {};
            for j = 1 : numOutliers
                currHighlight = cell2mat(fish(il).fishXClean(outlierIdx(j)));
                plot(TIME, currHighlight, 'color', 'r', 'LineWidth', 2);
                legendInfo{j}=['Trial # ', num2str(outlierIdx(j))]; % or whatever is appropriate
            end

            lgd = legend(legendInfo);
            lgd.Title.String = 'Outliers:';
            lgd.Location = 'southeast';
            lgd.NumColumns = 2;
            lgd.FontSize = 6;
            lgd.ItemTokenSize = [12, 6];

            % Display trial retention rate
            percentRemaining = round((numTrials - numOutliers) / numTrials * 100, 1);
            subtitle({['Z-Threshold: ', num2str(zThreshold)]; ...
                [' Max# of Outlier Peaks: ', num2str(maxAllowedOutlierPeaks)]; ...
                ['RMS-Threshold: ', num2str(zRmsThreshold)]; ...
                [num2str(percentRemaining), '% of Trials Remaining']});

            % Plot shuttle
            % sColor = [0, 1, 1, 0.9];
            % sLineWidth = 2;
            % plot(TIME, shuttle, 'color', sColor, 'LineWidth', sLineWidth);

            % Plot all the clean trials
            dataPoint = 'fishXClean';
            idxRange = 1: numTrials;
            lineWidth = 1.2;
            color = [0.2 0.2 0.6 0.4]; % dark blue
            plotSingleRepTimeDomain (fish, TIME, dataPoint, il, idxRange, plotName, yLim, color, lineWidth);

            % Save this plot for this luminance in the designated folder
            % Create folder to put figures
            % savePath = 'E:\Summer_2021\LIMBS_Lab\LIMBS_SU_2022\powerpoints\plots_for_presentations\TD_Clean_up_Threshold_Testing\11-15\';
            % thisFigFolder = [savePath, 'z_', zString, '_lim_', num2str(maxAllowedOutlierPeaks) '_zRMS_', ...
            %     zRmsString, '_acc_thresh_', zAccString, '_enum_', num2str(multiPeakThresh), '\'];
            % mkdir(thisFigFolder);
            saveas(gcf, [savePath, fishName, '_il_', num2str(il), '_TD.png']);
            saveas(gcf, [savePath, fishName, '_il_', num2str(il), '_TD.pdf']);

        end

        % Plot the average x-pos
        plotAvg = 0;
        if plotAvg == 1
            figure();
            hold on

            % Plot the shuttle
            TIME = 0:0.04:20-0.04;
            shuttlePoint = 'shuttleX';

            idxRange = 1;
            plotName = [fishName, ' Average x-pos with il = ', num2str(il)];
            yLim = [-0.04 0.04];

            sColor = [1, 0, 0, 0.4];
            sLineWidth = 2.5;

            plotSingleRepTimeDomain (fish, TIME, shuttlePoint, il, idxRange, plotName, yLim, sColor, sLineWidth);


            % Plot the trials
            dataPoint = 'xClean02Mean';
            idxRange = 1;
            lineWidth = 2;
            % color = colors{il}; % dark green

            plotSingleRepTimeDomain (fish, TIME, dataPoint, il, idxRange, plotName, yLim, [0 1 0.3], lineWidth);

            % Save this plot for this luminance in the designated folder
            saveas(gcf, [thisFigFolder, fishName, '_il_TD_mean_', num2str(il), '.png']);
        end

    end
    % Save the struct (not needed)
    % 
    % thisStructFolder = [savePath, 'z_', zString, '_lim_', num2str(maxAllowedOutlierPeaks) '_zRMS_', ...
    %     zRmsString, '_acc_thresh_', zAccString, '_enum_', num2str(multiPeakThresh), '\'];
    % 
    % % mkdir(thisStructFolder);
    % save([thisStructFolder, fishName, 'Struct'], 'fish');


end


function [mahalDistances, zScores, reValues, imValues] = calculateMahalZScores(fftValsAll)
[numTrials, numFreqs] = size(fftValsAll);
mahalDistances = zeros(numTrials, numFreqs);

reValues = zeros(numTrials, numFreqs);
imValues = zeros(numTrials, numFreqs);

for thisFreq = 1: numFreqs
    data = fftValsAll(:, thisFreq);
    re = real(data);
    im = imag(data);
    X = [re im];

    % [New] Use robustcov() directly
    [sig,mu,mah] = robustcov(X);

    % Calculate the sample mean and covariance
    mu = [mean(re) mean(im)];
    S = cov(re, im);

    % Calculate the Mahalanobis Distance to the sample distribution
    R = mvnrnd(mu, S, 1000);
    thisMahal = mahal(X, R);

    %[New] Swap out the mahal result with robustcov(X)
    %mahalDistances(:, thisFreq) = thisMahal;
    mahalDistances(:, thisFreq) = mah;
    reValues(:, thisFreq) = re;
    imValues(:, thisFreq) = im;

    % numTrials x 1 z-socre stats based on mahal distances
    zScores = zscore(mahalDistances);

end
end


function [fftValsAll, fftMagsAll] = calculateFFTbyTrials(inputCell)

% Set up
numTrials = size(inputCell, 2);
fftValsAll = zeros(numTrials, 12);
fftMagsAll = zeros(numTrials, 12);

for trNo = 1 : numTrials
    data = fixSmallTL(inputCell{trNo}, 10);

    % Parameters
    fs = 25;               % Sampling frequency (Hz)
    n = length(data);      % 500
    f = (0 : n - 1) * fs / n;  % Frequency axis
    
    % Compute FFT
    y = fft(data);
    yMag = abs(y) / (n/2) * 100;     % divide by (n/2), then convert from cm to m

    %% OLD 2022 code: this is not normalized correctly
    % yMag = abs(y);

    fs = 25;
    f = (0 : length(y) - 1) * fs / length(y);

    % Locate the correct peaks
    k = [2, 3, 5, 7, 11, 13, 19, 23, 29, 31, 37, 41];
    FREQ_AMPLITUDE_SCALAR = 0.05;
    freqData = k * FREQ_AMPLITUDE_SCALAR;

    fIdx = cast(freqData(:) / 0.05 + 1, "uint8");

    fftVals = y(fIdx);
    fftMags = yMag(fIdx);

    % Populate the big fft mag matrix (N x 12)
    fftValsAll(trNo, :) = fftVals;
    fftMagsAll(trNo, :) = fftMags;
end
end


% Helper: fixSmallTL
% Updated 07.20.2022 by Joy Yeh
% Fill up small tracking losses
%
% Params:
% data: the time-domain x value passed in. Might have small tracking loss
% window: the movmedian window size (10 is good)
% 
% Returns: 
% fixedData: filled gap with movmedian and applied a butterworth filter.

function [fixedData] = fixSmallTL(data, windowSize)
% fillmissing()
filled = fillmissing(data, 'movmedian', windowSize);

if sum(isnan(filled)) > 0
    fixedData = filled;
    return;
else
    %butterworth
    fs=25;
    fc = 3; % 10 Hz cutoff
    Wn = fc /(fs / 2); % Cut-off for discrete-time filter
    [b,a] = butter(2, Wn);
    fixedData = filtfilt(b,a, filled);
end
end

% Helper: calculateRmsZScores.m
% Updated 11/09/2022
% Pass in the x-pos data and shuttle.
% Returns the z-scores relative to other trials of the same condition (N x
% 1 array)
%
% Params:
% shuttle: the shuttle reference signal
% thisCell: the cell of x-pos data
%
% Returns:
% rms: the actual rmsError values of the trials (N x 1)
% zRms: the z-score values (N x 1)

function [rmsValues, zRms] = calculateRmsZScores(shuttle, thisCell)
numTrials = size(thisCell, 2);

% Populate the RMS field
rmsValues = [];
for p = 1 : numTrials
    
   
    % [OLD] the matlab rms() function seems to give NaN in some cases
    % rmsValues = [rmsValues; rms(shuttle - thisCell{p})];
    
    % [NEW] revert back to the definition of RMS();
    rmsValues = [rmsValues; nansum(abs((shuttle - thisCell{p}).^2))];
    
end

% numTrials x 1 z-socre stats based on rms
zRms = zscore(rmsValues);

end

% Helper: calculateFreqResponseValues.m
% Updated 11.09.2022 by Joy Yeh
% Calculate open- and closed-loop frequency responses and return arrays
%
% Params:
% shuttle: t-d shuttle input
% thisCell: the mean x-pos cell to work on
%
% Returns:
% GM, gain, phase, cp, cpGain, cpPhase: 1 x 12 complex/real double arrays
function [GM, gain, phase, cp, cpGain, cpPhase] = calculateFreqResponseValues(shuttle, thisCell)
k = [2, 3, 5, 7, 11, 13, 19, 23, 29, 31, 37, 41];

S = fftshift(fft(shuttle));
FM = fftshift(fft(cell2mat(thisCell)));
GM = zeros(size(k));
    
for j = 1:length(k)
    GM(1,j) =  FM((length(FM)/2)+1+k(j)) / S((length(FM)/2+1)+k(j));
end


gain = abs(GM);
phase = calibratePhase(rad2deg(unwrap(angle(GM))));

cp = GM ./ (1 - GM);
cpGain = abs(cp);
cpPhase = calibratePhase(rad2deg(unwrap(angle(cp))));

end

function [M] = calibratePhase(M)
    for i = 1 : length(M)
        if M(i) > 45
            M(i) = M(i) - 360;
        end

        if M(i) < -270
            M(i) = M(i) + 360;
        end
    end
   
end

% Helper: plotSingleRepTimeDomain()
% Updated 06.21.2022 by Joy Yeh
% plot the time-domain of a trial with the specified data, time scale, and
% plot title that indicates the current exp condition.
%
% Params:
% fish: the fish struct with all the data
% time: the time scale on the x-axis
% dataPoint: the object / point to plot.
% cond: the conductivity
% il: illumination
% idx: trial number
% plotName: the name of the plot describing what's happening.
% yLim: the y-limit of the axis
% color: plot color
%
% Returns: 
% A time-domain plot.

function plotSingleRepTimeDomain (fish, time, dataPoint, il, idxRange, plotName, yLim, color, lineWidth)
    hold on
    for j = idxRange
        data = fish(il).(dataPoint){j};
        plot(time, data, 'color', color, 'LineWidth', lineWidth);
        
    end
    
    xlim([0 20]);
    ylim(yLim);
    xlabel('time(s)');
    ylabel('fish position(m)');
    
    title(plotName);
end