close all; clear all;
%% Load Data Files
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
countN = 1; countP = 1;
for i=1:length(SignalFiles)
    load([Dir,'\',SignalFiles(i).name ])
    close all;
    if isfield(Data,'NaturalIntensity')
    NaturalPSTH(countN,:) = Data.Psth_sort{end}(1,1:50);
    countN = countN+1;
    else
    ProstheticPSTH(countP,:) = Data.Psth_sort{end}(1,1:50);
    countP = countP+1;
    end
end
MNPSTH = mean(NaturalPSTH,1); MPPSTH = mean(ProstheticPSTH,1);
NVis = size(NaturalPSTH,1); NNir = size(ProstheticPSTH,1);
Names = {['Natural(N=',num2str(NVis),')'],['Prosthetic(N=',num2str(NNir),')']};
%% Ploting
figure();
f1 = plot(linspace(-10,490,50),MNPSTH/max(MNPSTH),Color='g',LineWidth=2);
hold on
f2 = plot(linspace(-10,490,50),MPPSTH/max(MPPSTH),Color='r',LineWidth=2);
xlabel('Time [ms]'); ylabel('Normalized Spiking Rate');
xlim([-10 500]) ; ylim([0 1.01]);
axis square; box off;
set(gca,'color','none','FontSize',15)
legend([f1,f2],Names);
legend(EdgeColor='none',Color='none',Location='northeast'); 
%% FWHM
    dataVector = MNPSTH/max(MNPSTH);
    %dataVector = MPPSTH/max(MPPSTH);
    maxValue = max(dataVector); % Find the maximum value in the vector
    halfMax = maxValue / 2; % Calculate the half maximum value

    % Find the indices where the dataVector first and last reach half the maximum
    index1 = find(dataVector >= halfMax, 1, 'first');
    index2 = find(dataVector >= halfMax, 1, 'last');

    % Calculate the full width at half maximum
    fwhm = index2 - index1;

