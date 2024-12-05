% Tail plotting playground
% Does tail interpolation actually mess with the tail?
% What happened to Hope Il = 1, trial = 2, rep 1? The only valid trial?

close all
%% 1. Load the full body + rotated struct

abs_path = 'C:\Users\joy20\Folder\SP_2024\LIMBS Presentations\data\fish_structs_2024\';
out_path = 'C:\Users\joy20\Folder\SP_2024\LIMBS Presentations\Outputs\';
all_fish = load(fullfile(abs_path, 'data_clean_body.mat'), 'all_fish').all_fish;
% load([abs_path, 'data_clean_body.mat']);

fishNames = {'Hope', 'Len', 'Doris', 'Finn', 'Ruby'}; % consistent with SICB

p2m = 0.0004;
num_frames = 500;

Fs = 25;
Fc = 5; % cutoff frequencyin Hz
Wn = Fc/(Fs/2); % Cut-off frequency for discrete-time filter
[b, a] = butter(2, Wn); % butterworth filter parameters

% Define time domain data and parameters
T = 1/Fs; % Sampling period
L = 500; % Length of signal (500 frames)
t = (0:L-1)*T; % Time vector 1 x 500

%--------------------------------------------------------------------------

% range_fish = 1:5;
% range_il = 1;
% range_trial = 2;
% range_rep = 2;

% Create the FFT result structure
res = struct();
for i = 1 : 5
    h =  findobj('type','figure');
    n_fig = length(h);

    fish_name = fishNames{i};

    res(i).fish_name = fish_name;
    

    num_ils = numel(all_fish(i).luminance);

    res(i).luminances = struct();

    for il = 1: num_ils

        num_trials = numel(all_fish(i).luminance(il).data);
        count = 1;

        % These are for individual trials
        res(i).luminances(il).x_tail = {};
        res(i).luminances(il).y_tail = {};
        res(i).luminances(il).y_fft_amp = {};
        res(i).luminances(il).y_fft_vel = {};

        for trial_idx = 1 : num_trials %range_trial

            res(i).luminances(il).lux = all_fish(i).lux_values(il);

            f = all_fish(i).luminance(il).data(trial_idx);
            v = all_fish(i).luminance(il).data(trial_idx).valid_both;
            v_percentage = all_fish(i).luminance(il).data(trial_idx).valid_tail_percent;

            for rep = 1: 3 %range_rep

                % Grab the data
                valid = v(rep);

                if valid == 1
                    % Get meta data
                    res(i).luminances(il).trial_idx = trial_idx;
                    res(i).luminances(il).rep = rep;

                    v_percent = v_percentage(rep);
                    res(i).luminances(il).valid_percent = v_percent;

                    % Get the tail data
                    x_field = ['x_rot_rep', num2str(rep)];
                    y_field = ['y_rot_rep', num2str(rep)];
    
                    x_tail = all_fish(i).luminance(il).data(trial_idx).(x_field);
                    y_tail = all_fish(i).luminance(il).data(trial_idx).(y_field);
    
                    % only keep tail
                    x_tail = x_tail(:, 12);
                    y_tail = y_tail(:, 12);

                    % res(i).luminances(il).x_tail(count) = x_tail;
                    % res(i).luminances(il).y_tail(count) = y_tail;

                    % Append x_tail to the cell array
                    res(i).luminances(il).x_tail{end+1} = x_tail;
                    res(i).luminances(il).y_tail{end+1} = y_tail;


                    % 
                    % ------------------ GET FFT -------------------------

                   
                    % positions = y_tail;
                    % 
                    % positions = positions - mean(positions);
                    % positions = positions * p2m * 100; % unit in cm
                    % [f1, P1] = single_sided_spectra(positions, Fs);
                    % 
                    % window_size = 5; % Adjust as needed
                    % amp_smooth = movmean(P1, window_size);
                    % 
                    % delta_t = T; % Time difference between frames
                    % delta_position = diff(positions); % Differences in position
                    % velocity = delta_position / delta_t; % Velocity is change in position over time
                    % [f2, P2] = single_sided_spectra(velocity, Fs);
                    % vel_smooth = movmean(P2, window_size);
                    [f1, amp_smooth, f2, vel_smooth] = fftAmpAndVelocity(y_tail, Fs, T, p2m);
                    res(i).luminances(il).y_fft_amp{end+1} = amp_smooth;
                    res(i).luminances(il).y_fft_vel{end+1} = vel_smooth;
                    count = count + 1;
                end

                plot_3_figs = 0;
                if plot_3_figs == 1
                    % ------------------Plot FFT situations ------------------
                    main_figure = figure('Position', [100, 30, 500, 700]);
                    % set(main_figure, 'Visible', 'off');
                    subplot(311)
                    hold on

                    % with interpolation (saved in the struct)
                    % plot(1:num_frames, y_tail, 'Color', 'k');
                    % % no interpolation (just calculated it)
                    plot(1:num_frames, y_filt, 'Color', 'r');

                    % legend('clean y', 'with filt', 'Location', 'northeast');
                    xlabel('Frame');
                    ylabel('y-positions (pixel)'); % THERE ARE SOME SLIGHT DIFFERENCES!

                    title([num2str(v_percent), ' ', fish_name, 'Tail X FFT. Il = ', num2str(il), ', trial idx = ', ...
                        num2str(trial_idx), ', rep = ', num2str(rep), ', V = ', num2str(valid)]);

                    subplot(312)
                    plot(f1, amp_smooth, 'Color', 'b');
                    xlabel('Frequency (Hz)');
                    ylabel('Amplitude Smoothed (cm)');
                    xlim([0, 4]);
                    title('FFT Amplitude vs. Frequency');

                    subplot(313)
                    plot(f2, vel_smooth);
                    xlabel('Frequency (Hz)');
                    ylabel('Velocity Smoothed (cm/s)');
                    xlim([0, 4]);
                    title('FFT Velocity vs. Frequency');

                    % fig_out_path = [out_path, 'position_timelines\', fish_name,'\'];
                    fig_out_path = [out_path, 'tail_FFT\', fish_name,'\', num2str(il), '\'];
                    if ~exist(fig_out_path, 'dir')
                        mkdir(fig_out_path);
                    end

                    fig_out_filename = ['Tail_FFT_X_', fish_name, '_Il_', num2str(il), ...
                        '_trial_', num2str(trial_idx), '_rep_', num2str(rep), '_V_', num2str(valid), '_', num2str(v_percent), '.png'];

                    saveas(main_figure, [fig_out_path, fig_out_filename]);
                    disp(['SUCCESS: ', fig_out_filename, ' is saved.']);

                    % ---------------------------------------------------------
                end
            end

           
        end

        % Get the eman y-value and mean FFT for this il
        y_tail_mean = mean(cell2mat(res(i).luminances(il).y_tail), 2);
        
        res(i).luminances(il).y_tail_mean = y_tail_mean;

        [f1, amp_smooth_mean, f2, vel_smooth_mean] = fftAmpAndVelocity(y_tail_mean, Fs, T, p2m);
 

                
            

    end

end


% ---------- NOT GONNA RE-ROTATE

% Recover from x_data_raw (rotate again, NO interpolate)
re_rotate_positions = 0;
if re_rotate_positions == 1

    rotated_x = zeros(500, 12);
    rotated_y = zeros(500, 12);

    f = all_fish(fish_idx).luminance(il).data(trial_idx);
    x_field_raw = ['x_rep', num2str(rep)];
    y_field_raw = ['y_rep', num2str(rep)];

    for frame_idx = 1:500
        % Get the 12 points throughout the image
        x = f.(x_field_raw)(frame_idx, :);
        y = f.(y_field_raw)(frame_idx, :);

        % Linear fit with the first 3 points and rotate the fish
        coefficients = polyfit(x(1:3), y(1:3), 1);
        theta = atan(coefficients(1));
        [rotated_x(frame_idx, :), rotated_y(frame_idx, :)] = rotatePoints(x, y, x(2), y(2), theta);
    end

    % Only keep tail point for now
    rotated_x = rotated_x(:, 12);
    rotated_y = rotated_y(:, 12);
end


function [f1, amp_smooth, f2, vel_smooth] = fftAmpAndVelocity(X, Fs, T, p2m)
X = X - mean(X);
X = X * p2m * 100; % unit in cm
[f1, P1] = singleSidedSpectra(X, Fs);

window_size = 5; % Adjust as needed
amp_smooth = movmean(P1, window_size);

delta_t = T; % Time difference between frames
delta_position = diff(X); % Differences in position
velocity = delta_position / delta_t; % Velocity is change in position over time
[f2, P2] = singleSidedSpectra(velocity, Fs);
vel_smooth = movmean(P2, window_size);
end


% Helper: calculate FFT amplitude of input data X and the sampling freq Fs
function [f,P1] = singleSidedSpectra(X,Fs)

X(isnan(X))=[];
X = X - mean(X);

Y = fft(X);
L = length(X);

P2 = abs(Y/L * 2);
P1 = P2(1:round(L/2)+1);

f = Fs*(0:round(L/2))/L;
end


% Helper: rotate all the x and y data w.r.t the origin point
function [rotated_x, rotated_y] = rotatePoints(x, y, origin_x, origin_y, theta_rad)
% Rotate each point around the origin
rotated_x = (x - origin_x) * cos(-theta_rad) - (y - origin_y) * sin(-theta_rad) + origin_x;
rotated_y = (x - origin_x) * sin(-theta_rad) + (y - origin_y) * cos(-theta_rad) + origin_y;
end

