%% Fig_2ABCDEF.m
% Updated 09.20.2024
% LIMBS Lab
% Author: Huanying (Joy) Yeh

% Experiment Name: Eigenmannia Virescens Luminance + Locomotion Comparisons
%
% Content:
% - Generate plots for the body dynamics: "Rainbow" plots of points across
% the body, and the time-domain trajectories of these points
% - Display the following representative trials:
    % fish_idx = 1;
    % fish_name = 'hope';
    % il = 2, trial_number = 19, rep = 2 (0.4 lx)
    % il = 4, trial_number = 37, rep = 1 (2.0 lx)
    % il = 13, trial_number = 18, rep = 3 (150 lx)

% - Save most frames to "figures_archive"
%
% Output (for archive):
% "fig02a01_body_bending_vs_illuminance.png", with video frames overlaying the rainbow colors.
% "fig02a02_body_bending_timeline.png", it's the right panel with rainbow
% color timelines
% "fig02a03" and "fig02a04", same panels but for a different trial.

%% 1. Specify folder paths
parent_dir = fullfile(pwd);
abs_path = fullfile(parent_dir, 'data_structures\');

out_data_path = fullfile(parent_dir, 'data\media\');
if ~exist(out_data_path, 'dir')
    mkdir(out_data_path);
end

out_path = fullfile(parent_dir, 'figures\');
out_pdf_path = fullfile(parent_dir, 'figures_pdf\');

out_archive_path = fullfile(parent_dir, 'figures_archive\fig04a_body_dynamics\');
if ~exist(out_archive_path, 'dir')
    mkdir(out_archive_path);
end

pdf_path = fullfile(parent_dir, 'figures_pdf\');

close all

num_body_pts = 12;

%% 2. Load Video Frames

% Check validity
all_fish = load(fullfile(abs_path, 'data_clean_body.mat'), 'all_fish').all_fish;

fishNames = {'Hope', 'Len', 'Doris', 'Finn', 'Ruby'}; % consistent with SICB

p2m = 0.0004;
num_frames = 500;

%% 3. [FIRST TIME RUN] Copy the dedicated videos to the "data/media/" folder

%% Current Version 05/03
% Hope IL = 1, trial = 3, rep = 1 
% Hope IL = 9 Trial 6 Rep 3
% Hope Il = 14 Trial = 6 Rep 2

% frame_idx = 251 + rep * 500;
    % base_path = 'C:\Users\joy20\Folder\SP_2024\data\body_bending\';
    % dest_path01 = copy_videos(base_path, out_data_path, 'hope', 1, '03'); % rep 1, frame idx = 751
    % dest_path02 = copy_videos(base_path, out_data_path, 'hope', 4, '37'); % rep 3, frame idx = 1751
    % dest_path03 = copy_videos(base_path, out_data_path, 'hope', 13, '18');  % rep 2, frame idx = 1251


%% 4. Loop through and plot
for img_num = 1 : 3

    copy_videos_switch = 0; % already copied; turn this off

    if copy_videos_switch == 1
         %% 5. Save the video frame
    out_frame_filename = ['fig02a0', num2str(img_num), '_vid_frame_', fish_name, '_il_', num2str(il), ...
        '_trial_', num2str(trial_number), '_rep_', num2str(rep), '_frame_', num2str(frame_idx), '.png'];
    imwrite(frame, [out_path, out_frame_filename]);
    disp(['SUCCESS: ', out_frame_filename, ' is saved.']);

    end

    fish_idx = 1;
    fish_name = 'hope';

    % [NEW] 09/20/2024
     if img_num == 1
        il = 2;
        trial_number = 19;
        % dest_path = dest_path01;
        rep = 2;
        letter = 'a';

    elseif img_num == 2
        il = 4;
        trial_number = 37;
        % dest_path = dest_path02;
        rep = 1;
        letter = 'b';

    else
        il = 13;
        trial_number = 18;
        % dest_path = dest_path03;
        rep = 3;
        letter = 'c';
    end

    % Get video frames
    % frame_idx = 251 + rep * 500;
    % frame = display_vid_frame(dest_path, frame_idx); % TODO; replace this with original frame, not processed

   
    % Get the target data
    myColorMap = jet(num_body_pts);
    trial_idx = findFieldIdx(all_fish(fish_idx).luminance(il).data, trial_number);

    field_name_x = ['x_rot_rep', num2str(rep)];
    field_name_y = ['y_rot_rep', num2str(rep)];

    xData = all_fish(fish_idx).luminance(il).data(trial_idx).(field_name_x);
    yData = all_fish(fish_idx).luminance(il).data(trial_idx).(field_name_y);

    %% 6. Plot and save rainbow scatterplot
    title_text = ['Rainbow Scatterplot: ', fish_name, ' Il = ', num2str(il), ' Trial = ', num2str(trial_number), ' Rep = ', num2str(rep)];
    out_rainbow_fig = plotScatterWithColorbar(xData, yData, myColorMap, title_text);

    out_rainbow_filename = ['fig02a0', num2str(img_num), '_rainbow_', fish_name, '_il_', num2str(il), ...
        '_trial_', num2str(trial_number), '_rep_', num2str(rep), '.png'];
    saveas(out_rainbow_fig, [out_path, out_rainbow_filename]);

     out_rainbow_filename = ['fig02a0', num2str(img_num), '_rainbow_', fish_name, '_il_', num2str(il), ...
        '_trial_', num2str(trial_number), '_rep_', num2str(rep), '.pdf'];
    saveas(out_rainbow_fig, [out_pdf_path, out_rainbow_filename]);

    disp(['SUCCESS: ', out_rainbow_filename, ' is saved.']);


    %% 7. Plot and save timeline 
    time = 0:0.04:19.96;
    title_text = ['Timeline: ', fish_name, ' Il = ', num2str(il), ' Trial = ', num2str(trial_number), ' Rep = ', num2str(rep)];

    out_timeline_fig = plotTailPointTimeDomainScatterWithColorbar(time, yData, myColorMap, title_text);

    out_time_filename = ['fig02a0', num2str(img_num), '_timeline_', fish_name, '_il_', num2str(il), ...
        '_trial_', num2str(trial_number), '_rep_', num2str(rep), '.png'];
    saveas(out_timeline_fig, [out_path, out_time_filename]);

     out_time_filename = ['fig02a0', num2str(img_num), '_timeline_', fish_name, '_il_', num2str(il), ...
        '_trial_', num2str(trial_number), '_rep_', num2str(rep), '.pdf'];
    saveas(out_timeline_fig, [out_pdf_path, out_time_filename]);
    disp(['SUCCESS: ', out_time_filename, ' is saved.']);

end

%% Helper: copy videos to the data/media path
function destination_path = copy_videos(base_path, out_data_path, fish_name, il, trial_number)
path = [base_path, fish_name, '\', num2str(il), '\trial', trial_number, '*'];
match = dir(path);

source_path = fullfile(match.folder, match.name, "\vid_pre_processed.avi");
destination_path = [out_data_path, 'il_', num2str(il), fish_name, '_trial_', trial_number, '.avi'];
copyfile(source_path, destination_path);
disp("SUCCESS: Video is copied.");
end

%% Helper: display the given video frame from the input path
function frame = display_vid_frame(file_path, frame_idx)
vidReader = VideoReader(file_path);
frame = read(vidReader, frame_idx);
end

%% Helpers: given struct and the exp trial number, find its field index in the struct
function target_idx = findFieldIdx(input_struct, target_trial_idx)
leftmost_col = struct2cell(input_struct); %.trial_idx];
leftmost_col = cell2mat(squeeze(leftmost_col(1, :, :)));
target_idx = find(leftmost_col == target_trial_idx, 1);
end

%% Helper: Given x data, y data and color scheme, plot the scatter plot of all the body locations in the trial
function fig = plotScatterWithColorbar(xData, yData, myColorMap, title_text)
numPoints = size(xData, 2);

fig = figure('Position', [100, 100, 640, 190]);
set(fig, 'Color', 'none'); % Set background color to transparent

% set(gcf, 'Visible', 'off');

hold on;
for col = 1:numPoints
    scatter(xData(:, col), yData(:, col), 10, myColorMap(col, :), 'filled');
end
hold off;

colormap(myColorMap);

xticks([]);
yticks([]);

xlim([0, 640]);
ylim([-50, 250]);

title(title_text);
end

%% Helper: tail point t-d movement on the 20-second timeline
function fig = plotTailPointTimeDomainScatterWithColorbar(time, yData, myColorMap, titleText)

fig = figure('Position', [100, 100, 640, 190]);
% set(gcf, 'Visible', 'off');
numPoints = size(yData, 2); % Use all 12 points

% Plot all columns with the jet colormap
hold on;
for col = 1:numPoints
    plot(time, yData(:, col), 'Color', myColorMap(col, :), 'LineWidth', 2);
end
hold off;
colormap(myColorMap);

% Custom color bar
h = colorbar; 
custom_ticks = [0, 0.5, 1];
custom_labels = {'Head', 'FMiddle', 'Tail'};

title(h, 'Fish');
set(h, 'XTick', custom_ticks);
set(h, 'XTickLabel', custom_labels);

xlim([0, 20]); % timeline x-axis is out of 20 seconds
ylim([-50, 250]);

xlabel("Time(s)")
ylabel("Y-Pos (pixel)")

title(titleText);

grid on;
end


