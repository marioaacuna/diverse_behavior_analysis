% This script will detect and track the position of the animal in the
% 2-temperature arena for 2 animals. You need to select the compartment for
% each animal and then select the midline. Additionally, you can select if
% you want to see the video tracking or not. Finally, saves the data of percent
% time spent in each compartment in an excel file as well as trajectories and positions
% in a .mat file. 

cd('M:\Mario\Kristina')
% Some ideas taken from https://github.com/HanLab-OSU/MouseActivity/blob/master/README.md
%% Set up Vision system objects
clc
clear
close all


%% Select inputs
videos_rootpath = 'M:\Kristina\Conflicting motivation';
date_input = input('What date? YYMMDD with('') \n');
mat_filename = [videos_rootpath, '\', [date_input,'_Analysis_extracted_params.mat']];
xls_filename = [videos_rootpath, '\', [date_input,'_Analysis_timespent.xlsx']];


temp_name = questdlg('Select temperature', '','HOT', 'COLD', 'neutral', 'HOT');
switch temp_name
  case 'HOT'
    temp_path = [date_input,'-hot plate test 45deg'];
  case 'COLD'
    temp_path = [date_input,'-cold plate test 15deg'];
  case 'neutral'
    keyboard % to be added later maybe
end % switch


reply = questdlg('Do you want to diyplay the videos? Y/N');
% reply = input('Do you want to diyplay the videos? Y/N [Y]:','s');
if ~strcmp(reply, 'Yes')
    display_videos = 0;
else
    display_videos = 1;
end


videos_path = [videos_rootpath, '\', temp_path];
folder_components = dir(videos_path);
names = {folder_components.name};
movies = names(endsWith(names, '.mpg'));
movies = sort(movies)';

% Error in case no videos found
if isempty(movies)
    error('No videos were loaded. Check the name of the folder.')
end

ALL_DATA = struct();
%% Loop thu videos
all_data_cat =[];
DISTANCE_cat = [];
analysed_movies = cell(0,0);
ALL_positions = cell(0,0);
ALL_trajectories = cell(0,0);
is_roi_selected = 0;
for i_video = 1:length(movies)
    video_name = movies{i_video};
    total_path = [videos_path, '\' video_name];
    %% Read videos
%     vReader = vision.VideoFileReader(video_name, 'ImageColorSpace', 'YCbCr 4:2:2');
    the_vid = VideoReader(total_path); %#ok<TNMLP>
    vPlayer1 = vision.VideoPlayer('Position', [20,400,700,400], 'Name', [video_name, ': Animal TOP']);
    vPlayer2 = vision.VideoPlayer('Position', [800,400,700,400], 'Name', [video_name, ': Animal BOTTOM']);
    
    % Crop rectangle for each animal
    if ~is_roi_selected
        title_str = ('Draw ROI for each compartment. SELECT FIRST TOP ANIMAL');
        this_frame = read(the_vid, 100);
        ax = figure;
        hold on
        imshow(this_frame);
        title(title_str)
        hold off
        
        rect = NaN(2,4);
        for iid = 1:2
            rect(iid,:) = getrect(ax);
        end
        % Set midline
        this_frame_cropped =  imcrop(this_frame,rect(1,:));
        ax2 = figure;
        hold on
        imshow(this_frame_cropped)
        title ('Select now the midline that delimits each zone')
        hold off
        close(ax)
        % Select line
        [mid, midy] = getline(ax2);
        close(ax2)
        midline = mean(mid);
        is_roi_selected = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOT USED IN THE CURRENT FORM
    % % select Detector
    % detector_1 = vision.ForegroundDetector(...
    %     'NumGaussians', 3,...
    %     'NumTrainingFrames',50);%...
    % %     'MinimumBackgroundRatio',0.5);
    %
    % % select Detector
    % detector_2 = vision.ForegroundDetector(...
    %     'NumGaussians', 3,...
    %     'NumTrainingFrames',50);%...
    % %     'MinimumBackgroundRatio',0.5);
    %
    % min_area = 100;
    % blobAn_1 = vision.BlobAnalysis(...
    %     'BoundingBoxOutputPort', true,...
    %     'AreaOutputPort', 0,...
    %     'CentroidOutputPort', 0, ...
    %     'MinimumBlobArea', min_area);
    %
    % blobAn_2 = vision.BlobAnalysis(...
    %     'BoundingBoxOutputPort', true,...
    %     'AreaOutputPort', 0,...
    %     'CentroidOutputPort', 0, ...
    %     'MinimumBlobArea', min_area);
    % shapeInserter = vision.ShapeInserter('BorderColor','Black');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % set element structure
    str_1 = strel('square', 10);
    threshold = .20;
    objpixels = 2000; % min area for detecting a mouse
    %% Loop through frames
    tic
    frame_nr = 0;
    n_frames = the_vid.NumFrames;
    % centroid = cell(n_frames,1);
    % n_total_frames = the_vid.NumFrames;
    data_1 = NaN(0, 2); % 2 axes
    data_2 = data_1;
    disp('Processing videos, wait!')
    while hasFrame(the_vid) %frame_nr < n_frames
        frame_nr = frame_nr + 1;
        %     frame = vReader.step();
        frame = readFrame(the_vid);
        new_frame = imadjust(median(frame,3), [0.45 0.55],[]);
        new_frame = new_frame + .25;
        % select the animals
        animal_1 = imcrop(new_frame,rect(1,:)); %  [XMIN YMIN WIDTH HEIGHT];
        animal_2 = imcrop(new_frame, rect(2,:)); %  [XMIN YMIN WIDTH HEIGHT];
        
        new_frame_1 = bwareaopen(imbinarize(255-animal_1, threshold),objpixels) * 255; % remove any object smaller than 50 pixels
        new_frame_2 = bwareaopen(imbinarize(255-animal_2, threshold),objpixels) * 255; % remove any object smaller than 50 pixels
        mouse_1 = regionprops(new_frame_1,'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
        mouse_2 = regionprops(new_frame_2,'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
        
        %     areaArray = [mouse.Area];
        %     [~,idx] = max(areaArray);
        CentroidArray_1 = [mouse_1.Centroid];
        CentroidArray_2 = [mouse_2.Centroid];
        c1 = CentroidArray_1(~isnan(CentroidArray_1));
        c2 = CentroidArray_2(~isnan(CentroidArray_2));
        if isempty(c1) || isempty(c2)
            if frame_nr == 1 && ~isempty(c1) || ~isempty(c2)
                take_only_1_fov = true;
                fov_empty = find([isempty(c1), isempty(c2)]);
                frame_to_take = 2;
            else
                fov_empty = NaN;
                frame_to_take = frame_nr;
            end
            if isempty(c1)
                c1 = data_1(frame_to_take-1,:);
            else
                if isnan(fov_empty)
                    c2 = data_2(frame_to_take-1,:);
                else
                    c2 = [0 0];
                end
            end
        end
        data_1(frame_nr,:)= c1;
        data_2(frame_nr,:)= c2;
        
        %     mask_1 = detector_1.step(animal_1);
        %     mask_2 = detector_2.step(animal_2);
        %     new_mask_1 =  imopen(mask_1, str_1);
        %     new_mask_2 =  imopen(mask_2, str_1);
        %     bboxes_1 = blobAn_1.step(new_mask_1);
        %     bboxes_2 = blobAn_2.step(new_mask_2);
        out_1 = insertShape(animal_1,'FilledCircle',  [c1, 5], 'Color', 'green', 'LineWidth', 5);
        out_1_1 = insertShape(out_1, 'Line', [midline, 0, [midline 100]],'LineWidth', 2, 'Color', 'red');
        out_2 = insertShape(animal_2,'FilledCircle',  [c2, 5], 'Color', 'green', 'LineWidth', 5);
        out_2_2 = insertShape(out_2, 'Line', [midline, 0, [midline 100]],'LineWidth', 2, 'Color', 'red');
        
        %     mask = detector.step(new_frame);
        %     new_mask =  imopen(mask, str_1);
        % %     [~, centroid, bboxes] = blobAn.step(new_mask);
        %      bboxes = blobAn.step(new_mask);
        
        %     frame_markers = insertMarker(new_frame, centroid);
        %     out = insertShape(new_frame,'Rectangle',  bboxes, 'Color', 'green', 'LineWidth', 5);
        
        % Gather data
        %     if ~isempty(bboxes_1)
        %         data(frame_nr,1) = mean(bboxes_1(:,1));
        %     elseif ~isempty(bboxes_2)
        %         data(frame_nr,2) = mean(bboxes_2(:,1));
        %
        %     end
        if display_videos
            vPlayer1.step(out_1_1)
            vPlayer2.step(out_2_2)
        else
           waitbar(frame_nr / n_frames)
        end
%          f = waitbar(frame_nr / n_frames, 'Analysing this video, hold on');
    end
    %
%     release(vReader)
    release(vPlayer1)
    release(vPlayer2)
    toc
%     if ~display_videos
%         close(f)
%     end
    analysed_movies = [analysed_movies;video_name];
    disp('Done processing frames')
    
    %% analyse Time Spent in each compartment
    % Set midline
    
    propr_animal_1 = 100 *[sum(data_1(:,1) < midline) / frame_nr, sum(data_1(:,1) > midline) / frame_nr ];
    propr_animal_2 = 100* [sum(data_2(:,1) < midline) / frame_nr, sum(data_2(:,1) > midline) / frame_nr];
    data_left_1 = propr_animal_1(1);
    data_right_1 = propr_animal_1(2);    
    data_left_2 = propr_animal_2(1);
    data_right_2 = propr_animal_2(2);
    %% calculate trajectory
    FR = the_vid.FrameRate;
    total_time = sum(~isnan(data_1(:,1))) / FR; % seconds
%     x = sort(20*rand(100,1));
%     v = besselj(0,x);
    calibrate_pixels = 8.25 / diff(midy); % 1p = xcms
    time_F = linspace(0,1/total_time,sum(~isnan(data_1(:,1))));
    DISTANCE_all = [];
    trajectories = [];
    for iid  = 1:2
        switch iid
            case 1
                data_to_take = data_1(~isnan(data_1(:,1)),:);
                
            case 2
                data_to_take =  data_2(~isnan(data_1(:,1)),:);
        end
        
        xx = interp1(time_F,data_to_take(:,1), time_F,'spline');
        yy = interp1(time_F,data_to_take(:,2),time_F,'spline');
        differences = abs(diff(data_to_take(:,1))+diff(data_to_take(:,2)));
        distance = sum(differences); % check
        distance = distance * calibrate_pixels; % in cm 
        DISTANCE_all = [DISTANCE_all,  distance]; % check
        trajectories = [trajectories , [xx', yy']];
    end
    DISTANCE_cat = [DISTANCE_cat; DISTANCE_all]; % check
    positions = [data_1,data_2];
    %%
    data_cat = [data_left_1, data_right_1, data_left_2, data_right_2];
    all_data_cat = [all_data_cat ;data_cat];
    disp(['Animal TOP left/right: ',num2str(propr_animal_1 )])%, '%','; Right: ',num2str(100 - propr_animal_1), '%'])
    disp(['Animal BOTTOM left/right: ',num2str(propr_animal_2)])%, '%','; Right: ',num2str(100 - propr_animal_2), '%'])
    ALL_positions(i_video) = {positions};
    ALL_trajectories(i_video) = {trajectories};
    %% Terminate
    close all
    clear('vPlayer1')
    clear('vPlayer2')
    
end
DATA.(temp_name).positions = ALL_positions;
DATA.(temp_name).trajectories = ALL_trajectories;
DATA.(temp_name).movies = movies;
DATA.(temp_name).mazes = rect;
DATA.(temp_name).midline = [mid, midy]; 

if exist(mat_filename, 'file')
    DATA_pre  = load(mat_filename, 'DATA');
    prev_temp = fieldnames(DATA_pre.DATA);
    DATA.(cell2mat(prev_temp)) = DATA_pre.DATA.(cell2mat(prev_temp));
end

save(mat_filename, 'DATA')
% Write data to excel
data_table = array2table([all_data_cat, DISTANCE_cat] ,'VariableNames', {'Top left', 'Top right', 'Bottom left', 'Bottom right', 'Distance Top', 'Distance Bottom'}, 'RowNames', movies);
writetable(data_table, xls_filename, 'Sheet', temp_name, 'WriteRowNames', true)
disp(['all done for ', temp_name])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*AGROW>
