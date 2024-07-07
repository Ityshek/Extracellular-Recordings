clear all; clc; 
% Get a list of all files in the current directory
files = dir('C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Results\Data\MicronAlt');
%%
% Initialize variables to store dates and strings
dates = {}; Names = {};
Locations = {}; Params = {};
Channels = {};

% Loop through each file
for i = 1:length(files)
    filename = files(i).name;
    
    % Extract date and string from the file name
    % This depends on the format of your file names
    name = strsplit(filename,'_');
    date = name{end}(1:10); % replace with your extraction code
    Name = name{1}; Location = name{2}; % replace with your extraction code
    Channel = name{3}; Param = name{4};
    % Append the date and string to the lists
    dates{end+1} =  date; Names{end+1}  = Name;
    Locations{end+1} =  Location; Channels{end+1} = Channel;
    Params{end+1} =  Param;
end
%%
% Create a table from the dates and strings
T = table(dates', Names', Locations', Channels', Params');
%%
% Write the table to a new Excel file
writetable(T,...
    'C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Results\Data\AltUpdateList.xlsx');
