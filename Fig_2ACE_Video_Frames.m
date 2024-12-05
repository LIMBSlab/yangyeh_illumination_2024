%% Fig_2ACE_Video_Frames.m
% Updated 08.29.2024
% LIMBS Lab
% Author: Huanying (Joy) Yeh

% Experiment Name: Eigenmannia Virescens Luminance + Locomotion Comparisons
%
% Content:
% - Fig2_ACE: take snapshots of the video frames from "data/media/" folder

% Output (for archive):
% "fig06a01_body_bending_vs_illuminance.png", with video frames overlaying the rainbow colors.
% "fig06a02_body_bending_timeline.png", it's the right panel with rainbow
% color timelines
% "fig06a03" and "fig06a04", same panels but for a different trial.


close all

num_body_pts = 12;

%% 2. Load Video Frames
p2m = 0.0004;
num_frames = 500;

parent_dir = fullfile(pwd);
dest_path = fullfile(parent_dir, 'data_structures\video_frames\vid.avi');
% dest_path = "vid.avi";
% Get video frames

frame_indices = [1, 100, 200, 300, 400, 500, 1200, 1150];

for frame_idx = frame_indices
frame = display_vid_frame(dest_path, frame_idx); 

%% 5. [INPUT] Modify the file names, and then save the video frames
%% Overlaying happens in Affinity Designer - Check Fig_2ABCDEF.m for tracked datapoints
out_frame_filename = ['hope_il01_trial03_frame_', num2str(frame_idx), '.png'];
imwrite(frame, out_frame_filename);
disp(['SUCCESS: ', out_frame_filename, ' is saved.']);
end

%% Helper: display the given video frame from the input path
function frame = display_vid_frame(file_path, frame_idx)
vidReader = VideoReader(file_path);
frame = read(vidReader, frame_idx);
imshow(frame)
end

