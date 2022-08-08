function calculate_time_spent_center()
%% Preamble
% this script will take the results from the tracking script and will
% calculate trayectories and distance travelled in the open field.

clear 
clc
close all

do_figure = 1;

% length_maze_cm = 40;
% size_center
% length_center = 10; % 10 b 10cm
% rootpath = '\\pc150\groupnevian\Marta\Behavior\OpenField';
rootpath = 'M:\Federica\5HT-7\openfield';
% session = 'session 1 05.08.2021';
%session  = 'group 3';
xls_filename = [rootpath, '\', 'openfield.xlsx'];
results_path =  [rootpath,'\', 'results'];
file_path_figures = [rootpath, '\','Figures','\'];
% Read videos

all_parts = dir(results_path);
all_parts_names = {all_parts.name};
files_to_analyse = all_parts_names(endsWith(all_parts_names, '.mat'));

%%
n_files = length(files_to_analyse);
TIME_CENTER =[];
for ifile = 1:n_files
    this_file = files_to_analyse{ifile};
    video_str = strsplit(this_file, '.mat');
    video_str = video_str{1};
    % Load result parts 
    %     xls_filename = [results_path, '\', [video_str,'.xlsx']];
    mat_filename = [results_path, '\', [video_str,'.mat']];
    DATA = load(mat_filename, 'DATA');

    positions = DATA.DATA.positions{1, 1};
    
    minx = min(positions(:,1));
%     maxx = max(positions(:,1));
    
    miny = min(positions(:,2));
    maxy = max(positions(:,2));
    
    % bring the positions to 0 
    positions(:,1) = positions(:,1) - minx;
    positions(:,2) = positions(:,2) - miny;
    

    cal_line = DATA.DATA.cal_line;
    length_cal_line = cal_line(2) - cal_line(1) ; % in px
%     cal_px =  length_maze_cm / length_cal_line;
%     length_this_center = length_center / cal_px;
    
    % set square at 1/4 and 3/4 of the total length
%     adjusted_differece = maxy - length_cal_line;
    adjusted_differece  = 0;
    sc_diff = length_cal_line - 1/5* length_cal_line - adjusted_differece;
    sc_diff_1 =  1/5* length_cal_line - adjusted_differece;
    true_x = sc_diff_1 < positions(:,1) & positions(:,1) < sc_diff;
    true_y = sc_diff_1 < positions(:,2) & positions(:,2) < sc_diff;
    in_square_idx = (sum([true_x, true_y], 2) ==2);
    in_square = positions(in_square_idx, :);
    total_time_spent_in_center = 100 * length(in_square) /  length(positions); % in( (%)
    
    if do_figure
        %     test figure
        fig1 = figure;
        plot(positions(:,1), positions(:,2), 'o')
        hold on
        plot(in_square(:,1), in_square(:,2), 'ro')
        title(['percent : ', num2str(total_time_spent_in_center), ' - ', video_str], 'Interpreter', 'none')
        axis off
        if ~exist(file_path_figures, 'dir'), mkdir(file_path_figures), end
        filename_figure = [video_str,'.pdf'];
        saveas(fig1,[file_path_figures, filename_figure])
    end
    TIME_CENTER = [TIME_CENTER, total_time_spent_in_center]; %#ok<AGROW>
    
end

%% Write to excel
% 
% if startsWith(session, 'session 1')
%     toadd = '1';
% elseif startsWith(session, 'session 2')
%     toadd = '2';
% else
%     keyboard
% end
data_table = array2table(TIME_CENTER' ,'VariableNames', {'Percent time in center'}, 'RowNames', files_to_analyse');
writetable(data_table, xls_filename, 'Sheet', ['time_in_center'], 'WriteRowNames', true)

disp('all done')
end