clear all
close all
% Define path for data
Dir=uigetdir('*.fig','Select a Folder to Load Data Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = [Dir,'\'];

for i= 1:length(SignalFiles)
    load([pathname,SignalFiles(i).name])
    LatencyCut = 10*440; % Define the cut to split spike groups
    AA = 0;
    for r=1:length(Data.ind_rast_sort{end})
    AA = AA+size(Data.ind_rast_sort{end}{1,r},2);
    end


end