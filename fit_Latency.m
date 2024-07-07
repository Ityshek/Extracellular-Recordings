%% Load Pre-Proccesed Data
clear all
close all
% Define path for data
Dir=uigetdir('*.fig','Select a Folder to Load Data Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = [Dir,'\'];

% Create Variables
Intensity = []; StimDur = []; Latency = [];
NatTable = table(Intensity,StimDur,Latency);
ProsTable = table(Intensity,StimDur,Latency);

% Load Relevant Data into Tables
for i= 1:length(SignalFiles)

    load([pathname,SignalFiles(i).name])
    Intensity = []; StimDur = []; Latency = [];
    if isfield(Data,'NaturalLatency') % Load Data for Prosthetic Stimulation
    IsRespond =  find(~isnan(Data.NaturalLatency{1})); % Find the Intensity trials in which the Unit has responded
    Intensity = Data.NaturalIntensity(IsRespond);
    StimDur = ones(length(Data.NaturalIntensity(IsRespond)),1)*str2num(Data.StimDur);
    Latency = Data.NaturalLatency{1}(IsRespond);
    TempTable =   table(Intensity,StimDur,Latency);  
    NatTable = [NatTable ;TempTable];
    else % Load Data for Prosthetic Stimulation
    IsRespond =  find(~isnan(Data.ProstheticLatency{1})); % Find the Intensity trials in which the Unit has responded
    Intensity = Data.ProstheticIntensity(IsRespond);
    StimDur = ones(length(Data.ProstheticIntensity(IsRespond)),1)*str2num(Data.StimDur);
    Latency = Data.ProstheticLatency{1}(IsRespond);
    TempTable =   table(Intensity,StimDur,Latency);  
    ProsTable = [ProsTable ;TempTable];
    end
close all;
end

save([pathname 'ProccesedData\LatenctDataTables.mat'], 'NatTable', 'ProsTable')
%% Load Proccesed Data
clear all; close all;
[file,path]=uigetfile('*.mat','Select a Folder to Load Fig Files From');
load([path,file]);
%% Calculations
% Calc for Nat
AvgLat = [];
DursN = unique(NatTable.StimDur);
IntsN = unique(NatTable.Intensity); VarTypes = {'double','double','double'};
LatencyNat = table('Size',[11 3],'VariableTypes',VarTypes,'RowNames', string(IntsN),'VariableNames',string(DursN));
for i=1:length(DursN)
    for k=1:length(IntsN)
    Idxs = intersect(find(NatTable.StimDur == DursN(i)),find(NatTable.Intensity == IntsN(k)));
    AvgLat(k) = mean(NatTable.Latency(Idxs));
    end
LatencyNat(:,string(DursN(i))) = table(AvgLat');
end
% Calc for Pros
AvgLat = [];
DursP = unique(ProsTable.StimDur);
IntsP = unique(ProsTable.Intensity); VarTypes = {'double','double'};
LatencyPros = table('Size',[17 2],'VariableTypes',VarTypes,'RowNames', string(IntsP),'VariableNames',string(DursP));
for i=1:length(DursP)
    for k=1:length(IntsP)
    Idxs = intersect(find(ProsTable.StimDur == DursP(i)),find(ProsTable.Intensity == IntsP(k)));
    AvgLat(k) = mean(ProsTable.Latency(Idxs));
    end
LatencyPros(:,string(DursP(i))) = table(AvgLat');
end
%% Plotting
colors = ['r','g','b'];
figure();
subplot(1,2,1)
ax = gca;
for i=1:length(DursN)
    D = flip(DursN); % Sort Durs from Longest to Shortest
    Idxs = find(~isnan(LatencyNat.(string(D(i)))));
    FN = fit(str2double(LatencyNat.Row(Idxs)),LatencyNat.(string(D(i)))(Idxs),'poly1'); % Fit Latency Means without NaNs
    HN{i} = plot(FN,str2double(LatencyNat.Row(Idxs)),LatencyNat.(string(D(i)))(Idxs));
    set(HN{i},'color',colors(i),'DisplayName',[num2str(D(i)) 'ms Pulse'],'LineWidth',2,'MarkerSize',10)
    hold on
end
ax.FontSize = 15; %ax.XScale = "log";
xlabel('Intensity [nW/mm^2]'); ylabel('Latency [ms]')
ylim([0 200]); xlim([10 600])
legend(EdgeColor='none',Color='none');
axis square; box off; set(gca,'color','none')

title('Natural Vision')
subplot(1,2,2)
ax = gca;
for i=1:length(DursP)
    D = flip(DursP); % Sort Durs from Longest to Shortest
    Idxs = find(~isnan(LatencyPros.(string(D(i)))));
    FP = fit(str2double(LatencyPros.Row(Idxs))*1000,LatencyPros.(string(D(i)))(Idxs),'poly1'); % Fit Latency Means without NaNs
    HP{i} = plot(FP,str2double(LatencyPros.Row(Idxs))*1000,LatencyPros.(string(D(i)))(Idxs));
    set(HP{i},'color',colors(i),'DisplayName',[num2str(D(i)) 'ms Pulse'],'LineWidth',2,'MarkerSize',10)
    hold on
end
ax.FontSize = 15; %ax.XScale = "log";
xlabel('Intensity [\muW/mm^2]'); ylabel('Latency [ms]')
ylim([0 200]); xlim([10 3100])
legend(EdgeColor='none',Color='none'); 
axis square; box off; set(gca,'color','none')
title('Prosthetic Vision')