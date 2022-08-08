%% Preamble
% Plot heatmap of animal trajectory for a rectangular maze 
% cd('/Users/marioacuna/Desktop/') % change this for the actual directory
% clear
close all
clc
%%
cd('M:\Mario\Kristina\')
rootpath = 'M:\Kristina\Conflicting motivation';

indate = input('Select date: \n');
if ~ischar(indate)
    indate = num2str(indate);
end

%% Load data
inexp = input('Select experiment (COLD, HOT) in quatations (''): \n');

data = load([rootpath,'\',indate,'_Analysis_extracted_params.mat']);

%% get data per animal per video
this_data = data.DATA.(inexp).trajectories;
% Loop through videos
n_videos = size(this_data,2);
% Plot parametes
HSIZE = [10 10];SIGMA = 5;
grp_px = 125;
h = fspecial('gaussian', HSIZE,SIGMA);
names = data.DATA.(inexp).movies;
for iv = 1:n_videos
    % Loop through animals
    f_temp = figure('pos', [0 0 1800 1800], 'color', 'w');
    
    animal1 = this_data{1,iv}(:,[1,2]);
    try
        animal2 = this_data{1,iv}(:,[3,4]);
        n_animals = 2;
    catch
        n_animals = 1;
    end
    sp = [];
    for iid = 1:n_animals
        sp(iid) = subplot(2,1,iid);
        switch iid
            case 1
                this_animal_data = animal1;
            case 2
                this_animal_data = animal2;
        end
        %% get the 2d histogram
        N = hist2d(this_animal_data(:,1), this_animal_data(:,2), grp_px, 'countdensity');
        y = filter2(h, N);
        hi = imagesc(y);
%         xlim([0 500])
%         ylim([0 150])
        axis off 
        % load(KV)
        colormap(KV)% KV: variable created from colormap editor (run colormapeditor)
        title(['mouse ',num2str(iid)])
        %     set(hi,axis,off)
%         set(sp(iid), 'axis', 'square')
    end
    titlestr = [indate , ' - ',names{iv}, ' - ', inexp];
    sgtitle(titlestr)
    % save figure
    fig_path= [rootpath,'\', 'figures_passive' ,'\', inexp,'\' ];
    if ~exist(fig_path, 'dir')
        mkdir(fig_path)
    end
    fig_filename = [fig_path, titlestr, '.png' ];
    saveas(f_temp,fig_filename)
    %%
end