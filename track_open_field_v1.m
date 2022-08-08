% This script will detect and track the position of the animal in the
% Open Field for one animal. You need to select the compartment for
% each animal. Additionally, you can select if
% you want to see the video tracking or not. 

% Some ideas taken from https://github.com/HanLab-OSU/MouseActivity/blob/master/README.md
%% Set up Vision system objects
clc
clear
close all


%% Select inputs
rootpath = 'M:\Federica\5HT-7\openfield';
videos_rootpath = 'M:\Federica\5HT-7\openfield';
reply = questdlg('Do you want to display the videos? Y/N');
if ~strcmp(reply, 'Yes')
    display_videos = 0;
else
    display_videos = 1;
end

min_area = 450;
videos_path = [videos_rootpath];
folder_components = dir(videos_path);
names = {folder_components.name};
movies = names(endsWith(names, '.mp4'));
movies = sort(movies)';

ALL_DATA = struct();
%% Loop thu videos
analysed_movies = cell(0,0);
is_roi_selected = 0;
disp(movies)

video_name = input('What is the name of the video? ('') = ');
total_path = [videos_path, '\' video_name];
video_strs=  strsplit(video_name, '.');
video_str = video_strs{1};
xls_filename = [rootpath, '\results\', [video_str,'.xlsx']];
mat_filename = [rootpath, '\results\', [video_str,'.mat']];
figure_filename = [rootpath, '\results\', [video_str,'.pdf']];
time_start = input('In case your heands were on the video, write the time to start analysing here (in s), otherwise 0 -> ');
length_recording = input('What is the length of the recording (s)? -> ');

%% Read videos
the_vid = VideoReader(total_path); 
vPlayer1 = vision.VideoPlayer('Position', [20,400,700,400], 'Name', [video_name, ': Animal TOP']);

% Crop rectangle for each animal
if ~is_roi_selected
    title_str = ('Draw arena');
    this_frame = read(the_vid, 2000);
    ax = figure;
    hold on
    imshow(this_frame);
    title(title_str)
    hold off
    arena = getrect(ax);
    % Set midline
    this_frame_cropped =  imcrop(this_frame,arena);
    [l, w, ~] = size(this_frame_cropped);
    ax2 = figure;
    hold on
    imshow(this_frame_cropped)
    title ('Select now the center and calibration line')
    hold off
    close(ax)
    
    % Select center
    center= getrect(ax2);
    
    % Set calibration line
    [mid, midy] = getline(ax2);
    
    close(ax2)
end


% set element structure
str_1 = strel('square', 10);
threshold = .20;
objpixels = 2000; % min area for detecting a mouse
%% Loop through frames
tic
frame_nr = 1;
n_frames = the_vid.NumFrames;
frame_rate = round(the_vid.FrameRate);
from_frame = time_start * frame_rate;
frame_nr = frame_nr +  from_frame -1;
data_1 = NaN(length(from_frame-1 : n_frames), 2); % 2 axes
data_2 = data_1;
disp('Processing videos, wait!')
while frame_nr < n_frames % hasFrame(the_vid) %frame_nr < n_frames
    frame_nr = frame_nr + 1;

    frame = read(the_vid, frame_nr);
    
    frame = rgb2gray(frame);
    inputImage =   imcrop(frame, arena);
    [rows, columns] = size(inputImage);
    medianFilteredImage = medfilt2(inputImage);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         % automatic segmentation using kmeans
    %         numberOfClusters = 3;
    %         % Do kmeans clustering on the median filtered image.
    %         grayLevels = double(medianFilteredImage(:));
    %         [clusterIndexes, clusterCenters] = kmeans(grayLevels, numberOfClusters,...
    %             'distance', 'sqEuclidean', ...
    %             'Replicates', 1);
    %
    %         labeledImage = reshape(clusterIndexes, rows, columns);
    %         %         [maxValue, indexOfMaxValue] = max(clusterCenters);
    %         [minValue, indexOfMinValue] = min(clusterCenters);
    %
    %         mouse = labeledImage  == indexOfMinValue;
    %         mouse =  bwareafilt(mouse, 1);
    %         %         figure, imshow(mouse)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Traditional binarization and segmentation
    binary_image = imbinarize(255-medianFilteredImage,0.92);
    mouse = bwareaopen(binary_image, min_area) * 255; % remove any object smaller than mouse's area
    comp_1 = regionprops(mouse,'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
    CentroidArray_1 = [comp_1.Centroid];
    c1 = CentroidArray_1(~isnan(CentroidArray_1));
    this_area = [comp_1.Area];
    if isempty(c1) && isempty(this_area) % || isempty(c2)
        c1 = NaN;
        is_continue = 1;
    elseif this_area(end) <= min_area
        data_1(frame_nr,:)= NaN;
        is_continue = 1;
    else
        is_continue = 0;
    end
    
    % Output
    data_1(frame_nr,:)= c1;
    if is_continue
        out_1 = insertShape(medianFilteredImage,'FilledCircle',  [l/2,w/2, 10], 'Color', 'red', 'LineWidth', 5);
    else
        out_1 = insertShape(medianFilteredImage,'FilledCircle',  [c1, 5], 'Color', 'green', 'LineWidth', 5);
    end
    
    if display_videos
        vPlayer1.step(out_1)
    else
        waitbar(frame_nr / n_frames)
    end

end
%
release(vPlayer1)
toc
analysed_movies = [analysed_movies;video_name];
disp('Done processing frames')


%% calculate trajectory
% First, cut the recording to the length of the sesion (becuase the manual
% recording gives a bit more than 10 min for open field).

total_time = length(data_1) / frame_rate; % seconds
if total_time > length_recording % total length of recording in s
    frame_end = length_recording *frame_rate;
    data_diff = data_1(1:frame_end,:);
    total_time = length_recording;
else
    data_diff = data_1;
end
% remove nan values (is it necesary?)
% data_to_take = data_diff(~isnan(data_diff(:,1)),:);
data_to_take  = data_diff;
time_F = linspace(0,total_time,length(data_to_take));
xx = interp1(time_F,data_to_take(:,1), time_F,'linear');
yy = interp1(time_F,data_to_take(:,2),time_F,'linear');
%% calculate distance in pixels unsing euclidean distance between each point(frame)
distances = []; 
i_point = 1;
while i_point < length(data_to_take)
    % take point 0
    x0 = data_to_take(i_point,1);
    y0 = data_to_take(i_point,2);
    % the the next one
    x1 = data_to_take(i_point +1,1);
    y1 =  data_to_take(i_point+1 ,2);
    X = [x0,y0;x1,y1];
    d = pdist(X, 'euclidean');
    distances(i_point) = d;
    i_point = i_point + 1;
end
total_distance = nansum(distances);


%% Transform distance from px to cm
px = 40/diff(mid); % 40cm each side
calibrate_pixels = px; % 8.25 / diff(midy); % 1p = xcms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is without the euclidean distance calculated for each point (lines avobe)
% dist_E = sqrt(bsxfun(@minus,xx',xx).^2 + bsxfun(@minus,yy',yy).^2);
% distance = norm(xx'-yy');
% or
% differences = abs(diff(data_to_take(:,1))+diff(data_to_take(:,2))); % Do not take
% distance = sum(differences); % check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance = total_distance * calibrate_pixels; % in cm
trajectories = [xx', yy'];

positions = data_1;
%% Display distances
disp(['Animal moved : (cm) ',num2str(distance)])%, '%','; Right: ',num2str(100 - propr_animal_1), '%'])
ALL_positions = {positions}; % without removing nans
ALL_trajectories = {trajectories}; % interpolated to the total length of recording
%% Terminate
close all
clear('vPlayer1')


%% plot figure
final_figure = figure('Position', [600, 0, 600, 600]);
hold on
plot(w-xx, l-yy, 'ro', 'MarkerSize', 1);
box off, axis off
% suptitle('Trajectory of the animal in the Open field test')
title(['distance: ', num2str(distance)])
%%
% keyboard

DATA.positions = ALL_positions;
DATAtrajectories = ALL_trajectories;
DATA.movies = movies;
DATA.mazes = arena;
DATA.cal_line = mid; 

% if exist(mat_filename, 'file')
%     keyboard
%     DATA_pre  = load(mat_filename, 'DATA');
%     prev_temp = fieldnames(DATA_pre.DATA);
%     DATA.(cell2mat(prev_temp)) = DATA_pre.DATA.(cell2mat(prev_temp));
% end
%% save data

save(mat_filename, 'DATA')
saveas(final_figure,figure_filename)
% Write data to excel
data_table = array2table(distance ,'VariableNames', {'Distance (cm)'}, 'RowNames', {video_name});
writetable(data_table, xls_filename, 'WriteRowNames', true)
disp('all done')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*AGROW>
