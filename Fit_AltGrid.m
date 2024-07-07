% Fit Alternating Grid Results
close all; clear all;
%% Load Data Files
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = Dir;
Dates = strings; Channels = strings; 
Locs = strings; DurVec = strings;
StimCondition = []; Freqvec = [];
ResponseOrsVec = {}; ResponseCPMsVec = {};
OrsVec={}; CPMsVec = {};
for i=1:length(SignalFiles)
    load([pathname,'\',SignalFiles(i).name ])
    if regexp(SignalFiles(i).name,['Natural']) == 1 % Check for Natural / Prosthetic
        StimCondition(i) = 1;
    else
        StimCondition(i) = 2;
    end
    % Sort Trials by Unit
%     [a,b] = regexp(SignalFiles(i).name,'_[0123456789.]*');
%     SortVar(i).Date = SignalFiles(i).name(a+3:b-1);
    [c,d] = regexp(SignalFiles(i).name,'Loc[0123456789]*');
    SortVar(i).Loc = SignalFiles(i).name(d);
    [e,f] = regexp(SignalFiles(i).name,'Ch[0123456789]*');
    SortVar(i).Channel = SignalFiles(i).name(e+2:f);
    Channels = [Channels,string(SortVar(i).Channel)];
    Locs = [Locs,string(SortVar(i).Loc)];
    clear a b c d e f
    
    %SponMat = reshape(cell2mat(Data.Spon)',[size(Data.Response')])';
    %CombinedResponse = (Data.Response+Data.AltResponse)/2 - SponMat; % remove baseline activity
    %CombinedResponse = max(CombinedResponse,0); % turn any negative values to 0
    CombinedResponse = (Data.Response+Data.AltResponse)/2;
    
    % Find the maximum value in the matrix and its linear index
    [maxValue, linearIndex] = max(CombinedResponse(:));
    % Convert the linear index to row and column indices
    [MaxCPM, MaxOR] = ind2sub(size(CombinedResponse), linearIndex);
    
% Extract Response Vectors for Orientations and CPMs around maximal response
    ResponseOrsVec{i} = CombinedResponse(MaxCPM,:);
    ResponseCPMsVec{i} = CombinedResponse(:,MaxOR);
% Extract Orientations and CPMs vectors
OrsVec{i} = Data.Orientations; CPMsVec{i} = Data.CPM;
end
close all
%% Plot Individual Graphs
NVis = count(num2str(StimCondition),'1'); NNir = count(num2str(StimCondition),'2');
titles = {'Natural';'Prosthetic'}; LineStyles = {'-'; '-.'}; Names = {['Natural(N=',num2str(NVis),')'],['Prosthetic(N=',num2str(NNir),')']}; Colors = {'b','r'};
% Plot Orientation Response
OrientationsFig =  figure();
for i=1:length(ResponseOrsVec)
h{i} = plot(OrsVec{i},ResponseOrsVec{i}/max(ResponseOrsVec{i}),'LineStyle',LineStyles{StimCondition(i)},LineWidth=2,Marker='*',Color=Colors{StimCondition(i)})
hold on
end
xlabel('Orientation');ylabel('Normalized Response')
legend([h{1},h{end}],{'Natural';'Prosthetic'});
legend(EdgeColor='none',Color='none',Location='southeast'); 
axis square; box off;
set(gca,'color','none','FontSize',15)


% Plot CPM Response
CPMsFig = figure();
for i=1:length(ResponseOrsVec)
f{i} = plot(CPMsVec{i},ResponseCPMsVec{i}/max(ResponseCPMsVec{i}),'LineStyle',LineStyles{StimCondition(i)},LineWidth=2,Marker='*',Color=Colors{StimCondition(i)})
hold on
end
xlabel('CPM');ylabel('Normalized Response')
legend([f{1},f{end}],{'Natural';'Prosthetic'});
legend(EdgeColor='none',Color='none',Location='southeast'); 
axis square; box off;
set(gca,'color','none','FontSize',15)
%% Fit the Responses

% Define Gaussian function
gaussian = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2));
% Fit Orientations
OrsFitFig = figure();
for i=1:length(SignalFiles)
% Initial guess
b0 = [0.5; mean(OrsVec{i}); std(OrsVec{i})];
% Fit the Gaussian function to the data
%[beta{i},R{i},J{i},CovB{i},MSE{i}] = nlinfit(OrsVec{i}, ResponseOrsVec{i}/max(ResponseOrsVec{i}), gaussian, b0);
[beta{i},R{i},J{i},CovB{i},MSE{i}] = nlinfit(OrsVec{i},ResponseOrsVec{i}, gaussian, b0);

% Plot the fitted Gaussian curve
x_range = min(OrsVec{i}):5:179;
beta{i}(1) = 1;
%o{i} = plot(x_range, gaussian(beta{i}, x_range),'LineStyle','--',LineWidth=1,Color=Colors{StimCondition(i)});
o{i} = plot(x_range, gaussian(beta{i}, x_range),'LineStyle','-',LineWidth=2);%,Color=Colors{StimCondition(i)});
col = get(o{i},"Color");
hold on
%plot(OrsVec{i}, ResponseOrsVec{i}/max(ResponseOrsVec{i}),'LineStyle','none',Marker='*',MarkerSize=10,Color=Colors{StimCondition(i)});
plot(OrsVec{i}, ResponseOrsVec{i}/max(ResponseOrsVec{i}),'LineStyle','none',Marker='*',MarkerSize=10,Color=col);
hold on
end
xlabel('Orientation');ylabel('Normalized Response')
xlim([0 180]);
legend([o{1},o{end}],Names);
legend(EdgeColor='none',Color='none',Location='southeast'); 
axis square; box off;
set(gca,'color','none','FontSize',15)

% Fit CPMs
CPMsFitFig = figure();
for i=1:length(SignalFiles)
% Initial guess
Cb0 = [0.5; mean(CPMsVec{i}); std(CPMsVec{i})];
% Fit the Gaussian function to the data
%[Cbeta{i},CR{i},CJ{i},CCovB{i},CMSE{i}] = nlinfit(CPMsVec{i},(ResponseCPMsVec{i}/max(ResponseCPMsVec{i}))', gaussian, Cb0);
[Cbeta{i},CR{i},CJ{i},CCovB{i},CMSE{i}] = nlinfit(CPMsVec{i},(ResponseCPMsVec{i})', gaussian, Cb0);

% Plot the fitted Gaussian curve
x_range = 1:1:16.6;
Cbeta{i}(1) = 1;
%C{i} = plot(x_range, gaussian(Cbeta{i}, x_range),'LineStyle','-',LineWidth=1.5,Color=Colors{StimCondition(i)});
C{i} = plot(x_range, gaussian(Cbeta{i}, x_range),'LineStyle','-',LineWidth=2);%,Color=Colors{StimCondition(i)});
col = get(C{i},"Color");
hold on
%plot(CPMsVec{i}, ResponseCPMsVec{i}/max(ResponseCPMsVec{i}),'LineStyle','none',Marker='*',MarkerSize=12,Color=Colors{StimCondition(i)});
plot(CPMsVec{i}, ResponseCPMsVec{i}/max(ResponseCPMsVec{i}),'LineStyle','none',Marker='*',MarkerSize=12,Color=col);
hold on
end
xlim([0 15]); ylim([0 1.01]);
xlabel('CPM');ylabel('Normalized Response')
legend([C{1},C{end}],Names);
legend(EdgeColor='none',Color='none',Location='northeast'); 
axis square; box off;
set(gca,'color','none','FontSize',15)



