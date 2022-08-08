% This script will detect and track the position of the animal in the
% place aversion chamber for one animal. You need to select the compartment for
% each animal. Additionally, you can select if
% you want to see the video tracking or not. Finally, saves the data of percent
% time spent in each compartment in an excel file as well as trajectories and positions
% in a .mat file.

% cd('M:\Mario\Liselot')
% Some ideas taken from https://github.com/HanLab-OSU/MouseActivity/blob/master/README.md
%% Set up Vision system objects
clc
clear
close all


%% Select inputs
type_experiment = input('type experiment (CPA/GsDREADD)? (('')) = ');
date_experiment = input('date of experiment? ("yy.mm.dd", with quotation marks('')) = ');


% date_experiment= '21.07.28';
switch type_experiment
    case 'CPA'
    videos_rootpath = ['M:\Liselot\Data\7. Conditioned Place Avoidance Paradigm\', date_experiment];
    min_area = 7500;
    length_shortest_side = 19;

    case 'GsDREADD'
    videos_rootpath = ['M:\Liselot\Data\2. opto-5HT7\CPP\GsDREADD_CCI\', date_experiment];
    min_area = 500;
    length_shortest_side = 29.3;
end
reply = questdlg('Do you want to display the videos? Y/N');
% reply = input('Do you want to display the videos? Y/N [Y]:','s');
if ~strcmp(reply, 'Yes')
    display_videos = 0;
else
    display_videos = 1;
end

% min_area = 7500;
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
% for i_video = 1:length(movies)
%     video_name = movies{i_video};
total_path = [videos_path, '\' video_name];
xls_filename = [videos_rootpath, '\', [video_name,'.xlsx']];
mat_filename = [videos_rootpath, '\', [video_name,'.mat']];

time_start = input('If you know the time when you open the gates, write it here (in s) -> ');

%% Read videos
%     vReader = vision.VideoFileReader(video_name, 'ImageColorSpace', 'YCbCr 4:2:2');
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
    ax2 = figure;
    hold on
    imshow(this_frame_cropped)
    title ('Select now the midline that delimits each zone')
    hold off
    close(ax)
    
    % Select line
    [mid, midy] = getline(ax2);
  
    [l, w, ~] = size(this_frame_cropped);
    is_vertical = l > w;  % is vertical?
    
    if is_vertical
        midline = mean(midy);
        cal_line = mid;
    else
        midline = mean(mid);
        cal_line = midy;
    end
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

    %     frame = vReader.step();
    %         frame = readFrame(the_vid);
    frame = read(the_vid, frame_nr);
    
    frame = rgb2gray(frame);
    inputImage =   imcrop(frame, arena);
    [rows, columns] = size(inputImage);
    medianFilteredImage = medfilt2(inputImage);
    %         medianFilteredImage = 1 - medianFilteredImage;
    
    
    
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
    binary_image = imbinarize(255-medianFilteredImage,0.75); %0.92
    mouse = bwareaopen(binary_image, min_area) * 255; % remove any object smaller than mouse's area
    %         new_frame_1 = bwareaopen(mouse);
    %         new_frame_1 = bwareaopen(imbinarize(255-mouse, threshold),objpixels) * 255; % remove any object smaller than 50 pixels
    %         new_frame_2 = bwareaopen(imbinarize(255-animal_2, threshold),objpixels) * 255; % remove any object smaller than 50 pixels
    comp_1 = regionprops(mouse,'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
    %         comp_2 = regionprops(new_frame_2,'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
    
    %     areaArray = [mouse.Area];
    %     [~,idx] = max(areaArray);
    CentroidArray_1 = [comp_1.Centroid];
    %         CentroidArray_2 = [mouse_2.Centroid];
    c1 = CentroidArray_1(~isnan(CentroidArray_1));
    %         c2 = CentroidArray_2(~isnan(CentroidArray_2));
    this_area = [comp_1.Area];
    if isempty(c1) && isempty(this_area) % || isempty(c2)
        c1 = NaN;
        is_continue = 1;
        %                 c2 = data_2(frame_nr-1,:);
    elseif this_area(end) <= min_area
        data_1(frame_nr,:)= NaN;
        is_continue = 1;
    else
        is_continue = 0;
    end
    
    % Output
    data_1(frame_nr,:)= c1;
    if is_continue
        out_1 = insertShape(medianFilteredImage,'FilledCircle',  [w/2, l/2, 10], 'Color', 'red', 'LineWidth', 5);
    else
        out_1 = insertShape(medianFilteredImage,'FilledCircle',  [c1, 5], 'Color', 'green', 'LineWidth', 5);
    %         out_1_1 = insertShape(out_1, 'Line', [midline, 0, [midline 100]],'LineWidth', 2, 'Color', 'red');
    %         out_2 = insertShape(comp_2,'FilledCircle',  [c2, 5], 'Color', 'green', 'LineWidth', 5);
    %         out_2_2 = insertShape(out_2, 'Line', [midline, 0, [midline 100]],'LineWidth', 2, 'Color', 'red');
    end
    if display_videos
        vPlayer1.step(out_1)
        %             vPlayer2.step(out_2_2)
    else
        waitbar(frame_nr / n_frames)
    end
    

    %          f = waitbar(frame_nr / n_frames, 'Analysing this video, hold on');
end
%
release(vPlayer1)
%     release(vPlayer2)
toc
analysed_movies = [analysed_movies;video_name];
disp('Done processing frames')

%% analyse Time Spent in each compartment
% Set midline
if is_vertical
    data_axis = data_1(~isnan(data_1(:,1)),2);
else
    data_axis = data_1(~isnan(data_1(:,1)),1);
end

frame_to_use = frame_nr - from_frame;
frames_to_normalize = sum(~isnan(data_1(:,1)));
propr_animal_1 = 100 *[sum(data_axis < midline) / frames_to_normalize, sum(data_axis > midline) / frames_to_normalize ]; % [up, down] /  [left,right]
%     propr_animal_2 = 100* [sum(data_2(:,1) < midline) / frame_nr, sum(data_2(:,1) > midline) / frame_nr];
data_left_1 = propr_animal_1(1);
data_right_1 = propr_animal_1(2);
% data_left_2 = propr_animal_2(1);
% data_right_2 = propr_animal_2(2);

disp(['proportion L,R or U,D = ', num2str(propr_animal_1)])
%% calculate trajectory
FR = the_vid.FrameRate;
total_time = sum(~isnan(data_1(:,1))) / FR; % seconds
%     x = sort(20*rand(100,1));
%     v = besselj(0,x);
calibrate_pixels = length_shortest_side / diff(cal_line); % 1p = xcms




time_F = linspace(0,1/total_time,sum(~isnan(data_1(:,1))));
data_to_take = data_1(~isnan(data_1(:,1)),:);
xx = interp1(time_F,data_to_take(:,1), time_F,'spline');
yy = interp1(time_F,data_to_take(:,2),time_F,'spline');
%% Calculate distance moved by sum of single point distances
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

% differences = abs(diff(data_to_take(:,1))+diff(data_to_take(:,2)));
% distance = sum(differences); % check
distance = total_distance * calibrate_pixels; % in cm
DISTANCE_all = total_distance; % check
trajectories = [xx', yy'];

DISTANCE_cat = DISTANCE_all; % check
positions = [data_1];
%%
data_cat = [data_left_1, data_right_1];
all_data_cat =data_cat;
disp(['Animal left/right (or up/down): ',num2str(propr_animal_1 )])%, '%','; Right: ',num2str(100 - propr_animal_1), '%'])
% disp(['Animal BOTTOM left/right: ',num2str(propr_animal_2)])%, '%','; Right: ',num2str(100 - propr_animal_2), '%'])
ALL_positions = {positions};
ALL_trajectories = {trajectories};
%% Terminate
close all
clear('vPlayer1')
% clear('vPlayer2')


%% plot figure
final_figure = figure('Position', [600, 0, w, l]);
hold on
plot(w-data_1(:,1), l-data_1(:,2), 'ro', 'MarkerSize', 1);
%%

DATA.positions = ALL_positions;
DATAtrajectories = ALL_trajectories;
DATA.movies = movies;
DATA.mazes = arena;
DATA.midline = [mid, midy]; 

if exist(mat_filename, 'file')
    keyboard
    DATA_pre  = load(mat_filename, 'DATA');
    prev_temp = fieldnames(DATA_pre.DATA);
    DATA.(cell2mat(prev_temp)) = DATA_pre.DATA.(cell2mat(prev_temp));
end

save(mat_filename, 'DATA')
% Write data to excel
data_table = array2table([all_data_cat, DISTANCE_cat] ,'VariableNames', {'Top / left', 'Bottom / right', 'Distance (cm)'}, 'RowNames', {video_name});
writetable(data_table, xls_filename, 'WriteRowNames', true)
disp('all done')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*AGROW>
