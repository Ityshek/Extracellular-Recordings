function [Im, ImageFiles,Dir]=load_image_directory 
Dir=uigetdir('*.bmp','Select Folder to Load stimulus Images From');
% Dirmeta=uigetdir(Dir,'*.xml');
ImageFiles=dir(fullfile(Dir, '*.bmp'));
% meta_data=dir(fullfile(Dirmeta, '*.xml'));

for k=1:length(ImageFiles)
        Im{k}=(imread([Dir '\' ImageFiles(k).name]));
end