% Fig_S5.m
% Script: Time_Domain_Clean_Up_Mahal_RMS.m
% Add validity tags to fish structs
% Use Complex Mahalanobis Z-scores and TD x-pos RMS Error Z-scores.
% Use pre-cleaned data structures

% Updated 12.05.2024

%% 1. Specify folder paths
parent_dir = fullfile(pwd);
abs_path = fullfile(parent_dir, 'data_structures\');

out_path = fullfile(parent_dir, 'figures\');
out_pdf_path = fullfile(parent_dir, 'figures_pdf\');

out_archive_path = fullfile(parent_dir, 'figures_archive\fig04a_body_dynamics\');
if ~exist(out_archive_path, 'dir')
    mkdir(out_archive_path);
end

pdf_path = fullfile(parent_dir, 'figures_pdf\');

close all

fishNameList = {'hope', 'len', 'doris', 'finn', 'ruby'};
N_FISH = 5;

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
for idx = 4
    fishName = fishNameList{idx};

    % [NEW] the following line replaced "BigStruct"
    % load(['../Fish_Structs/', fishName, 'BigStruct.mat']);
    % load([pwd, '\11-15_Mahal_RMS_Tuned\Original\', fishName, 'Original.mat']);
    load([abs_path, '\pre_cleaned_original_structs\', fishName, 'Original.mat']);
    load([abs_path, '\pre_cleaned_original_structs\shuttle.mat']);
    % load('../shuttle.mat');

    numConditions = size(fish, 2);

    for il = 1 % 1 : numConditions

        thisCell = fish(il).fishXClean;
        numTrials = size(thisCell, 2);

        [fftValsAll, fftMagsAll] = calculateFFTbyTrials(thisCell);
        % 15x12, 15x12
        [mahalDistances, Zmahal, re, im] = calculateMahalZScoresFinalPlotting(fftValsAll);


        % Plot for one frequency on the Mahalanobis complex plane
        for freq_band = 2
            % Define the frequency band target
            Zmahal_target = Zmahal(:, freq_band);

            % Create a new figure
            fig = figure();
            set(fig, 'defaultLegendAutoUpdate', 'off');
            hold on

Zmahal_target = Zmahal(:, freq_band);

% Create a new figure
set(fig, 'defaultLegendAutoUpdate', 'off');
hold on

% Generate Gaussian distributed points using the sample mean and covariance
numGaussianPoints = 1000;  % Number of additional Gaussian points to plot (for smoother contours)
reMean = mean(re(:, freq_band));
imMean = mean(im(:, freq_band));
S = cov(re(:, freq_band), im(:, freq_band));
R = mvnrnd([reMean, imMean], S, numGaussianPoints);

% Calculate Mahalanobis distances for each point
mahalDistances = mahal(R, R);  % Mahalanobis distances of generated points

% Generate a grid for contour plotting
xrange = linspace(min(R(:,1)), max(R(:,1)), 100);
yrange = linspace(min(R(:,2)), max(R(:,2)), 100);
[X, Y] = meshgrid(xrange, yrange);

% Calculate the Mahalanobis distance for each point on the grid
invS = inv(S);  % Inverse of covariance matrix
Z_grid = arrayfun(@(x, y) sqrt(([x y] - [reMean imMean]) * invS * ([x y] - [reMean imMean])'), X, Y);  % Mahalanobis distance for each grid point

% Define the Z-score levels for contours
z_levels = [1, 1.5, 2.5];

% Plot the contours at specified Z levels and fill with magenta color
contour(X, Y, Z_grid, z_levels, 'LineColor', 'm');  % Draw contour lines at Z = 1, 1.5, and 2.5
for i = 1:length(z_levels)
    % Get the contour at the specific Z level
    [C, h] = contour(X, Y, Z_grid, [z_levels(i) z_levels(i)], 'LineStyle', 'none');
    
    % Extract the contour data from C matrix
    k = 1;
    while k < size(C, 2)
        num_points = C(2, k);  % Number of points in this contour
        h1 = fill(C(1, k+1:k+num_points), C(2, k+1:k+num_points), [1 0 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % Fill contour with magenta color (30% opacity)
        k = k + num_points + 1;
    end
end

% Now, plot the actual data points on top of the level sets
for i = 1:size(Zmahal_target,1)
    if Zmahal_target(i) > 2.5
        % Plot outliers (z > 2.5) in red
        h2 = plot(re(i, freq_band), im(i, freq_band), 'ro', 'MarkerFaceColor', 'r');
        % disp(['plotting trial #', num2str(i), ' in red, distance is ', num2str(mahalDistances(i, freq_band))])
    else
        % Plot normal points (z <= 2.5) in blue
        h3 = plot(re(i, freq_band), im(i, freq_band), 'bo', 'MarkerFaceColor', [0.2 0.2 0.6]);
    end
end

hold off;

% Label the axes
xlabel('Real Part');
ylabel('Imaginary Part');
title([fishName, ' il = ', num2str(il), ', Mahalanobis Distances, freq band = ', num2str(freq_band)]);

% Add a legend with specific labels for each plot
if exist('h2', 'var')
    legend([h1, h2, h3], 'level sets (Z = 1, 1.5, 2.5)', 'outlier datapoints', 'datapoints', 'Location', 'best');
else
    legend([h1, h3], 'level sets (Z = 1, 1.5, 2.5)', 'datapoints (z <= 2.5)', 'Location', 'best');
end

            % Prepare for saving the plot
            
            saveas(gcf, [out_path, fishName, '_il_', num2str(il), '_freq_', num2str(freq_band), '_MAHAL.png']);
            saveas(gcf, [out_pdf_path, fishName, '_il_', num2str(il), '_freq_', num2str(freq_band), '_MAHAL.pdf']);
        end

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

    % % Prepare for plotting
    % out_path = 'C:\Users\joy20\Folder\SP_2024\LIMBS Presentations\FA_2024\';
    % saveas(gcf, [thisFigFolder, fishName, '_il_', num2str(il), '_MAHAL.png']);

end

function [mahalDistances, zScores, reValues, imValues] = calculateMahalZScoresFinalPlotting(fftValsAll)
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
    y = fft(data);
    
    %% OLD 2022 code: this is not normalized correctly
    % yMag = abs(y);

    yMag = abs(y) / size(data, 1);

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

% Helper: calibratePhase
% Updated 08.23.2022 by Joy Yeh
% Add or subtract the phase by multiples of 2 pi (360 degrees).
%
% Params:
% M: the original phase matrix in degrees
% G: the complex matrix
%
% Returns:
% M: the phase matrix with each of the possible 12 phases mofified, if
% needed.

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

