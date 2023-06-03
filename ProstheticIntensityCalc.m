function [Data] = ProstheticIntensityCalc(Data,PSTHbinsize,fname,window,PSTH)
%% Prosthetic Intensity Response Curve
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if answer{1} == '1'
    Data.ProstheticIntensity = []; Data.ProstheticIntensityResponse = cell(1,length(Data.Clusters)); Data.ProstheticIntensityCount = cell(1,length(Data.Clusters));   Data.ProstheticLatency = cell(1,length(Data.Clusters));
    Data.StimThresh = cell(1,length(Data.Clusters));
end

% Calculation
% Amplitude to Intensity Conversion
B = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7;...
    0.42*10^(-3),0.3,0.9,1.45,1.93,2.38,2.8,3.15,3.6,4.05,4.47,4.88,5.3,5.83]; % 1st row are current in Amp. Second row are intensity in mW.
D = [0.8,0.95,1.25,1.6,2.1,2.7,3.9,4.9;...
    0.097,0.27,0.53,1.01,1.47,2.09,3.03,4.05];
%a = [strfind(fname,'Hz')+1,strfind(fname,'amp')];
%Amp = fname(a(1)+2:a(2)-1);
%Data.ProstheticIntensity = [Data.ProstheticIntensity;D(find(D == str2num(Amp))+1)]; %Convert to Intensity/mm^2 for 40% duty cycle.
% Data.ProstheticIntensity = [Data.ProstheticIntensity;B(find(B == str2num(Amp))+1)]; %Convert to Intensity/mm^2 for 40% duty cycle.
[a,b] = regexp(fname,'_[.0123456789]*Int');
Data.ProstheticIntensity = [Data.ProstheticIntensity;str2num(fname(a+1:b-3))]; % Stimulus Intensity in mW/mm^2.


% Response Calculation
for i=1:length(Data.Clusters)
    ResponseWindow = [(window(1)/10)+2:(window(2)/10)+2]; % Define time window for Prosthetic response (10-210ms post trigger). add 2 bins for -10 and 0 bins in PSTH.
    Data.ProstheticIntensityResponse{i} = [Data.ProstheticIntensityResponse{i};round(max(PSTH(ResponseWindow)),2)];
    Data.ProstheticIntensityCount{i} = [Data.ProstheticIntensityCount{i};Data.SpikeCount{i}];
end

%% Stimulation Threshold
for i=1:length(Data.Clusters)
    %[RS,RM] = std(Data.Psth_sort{Data.Clusters(i)},'omitnan'); % Calc std and mean of total acticvity
    [Val,Ind]= max(PSTH((window(1)/10)+3:(window(2)/10)+2));
    %R = mean(Data.Psth_sort{Data.Clusters(i)}(Ind+1:Ind+4));
    R = Data.SpikeCount{Data.Clusters(i)} ;
    RM = Data.Spon{end}; RS = Data.SponStd{end};
    if length(find(PSTH(window(1)/10:window(2)/10)>RM+1.5*RS)) > 3 % Check for at least 3 bins in PSTH higher then 95% of baseline activity inside response window.
        % if R > RM+2*RS % Look if R window is bigger then 95% of total activity
        Data.StimThresh{i} = [Data.StimThresh{i};1];

        % Calc Latency (Only if responsive)
        %         Latency = (find(PSTH((window(1)/10)+2:(window(2)/10)+2)...
        %             == max(PSTH((window(1)/10)+2:(window(2)/10)+2)),1))*PSTHbinsize;
        %         Data.ProstheticLatency{i} = [Data.ProstheticLatency{i};Latency];
        Latency = (min(find(PSTH((window(1)/10)+3:(window(2)/10)+2)>0.5*Val))+1)*10; % Calc Latency as first bin in smoothed PSTH to pass 50% from max response in window.
        Data.ProstheticLatency{i} = [Data.ProstheticLatency{i};Latency];
    else
        Data.StimThresh{i} = [Data.StimThresh{i};0];
        Data.ProstheticLatency{i} = [Data.ProstheticLatency{i};nan];
    end
end
% Prosthetic Response Latency
% for i=1:length(Data.Clusters)
%     V = min(find(round(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow),2) == Data.ProstheticIntensityResponse{i}(end))+min(ResponseWindow)-1);
%     Data.ProstheticLatency{i} = [Data.ProstheticLatency{i};(V-2)*PSTHbinsize]; % Calculate latency of response in mSec.
% end
if answer{2} =='1'
    % Plotting
    %   Spiking Rate
    for i=1:length(Data.Clusters)
        Data.IntensityCurve{i} = figure();

        %         spon = ones(1,length(Data.ProstheticIntensityResponse{i}))*(mean(Data.Spon{i})+(3*std(Data.Spon{i})));
        %         plot([Data.ProstheticIntensity],[Data.ProstheticIntensityResponse{i}],'-',Data.ProstheticIntensity,spon,'--')
        semilogx([Data.ProstheticIntensity],[Data.ProstheticIntensityCount{i}],'b',...
            [Data.ProstheticIntensity],ones(1,length(Data.ProstheticIntensityCount{i}))*mean(Data.Spon{i}),'r')
        xticks(linspace(0,round(max(Data.ProstheticIntensity),1),5))
        ylabel('Spike Count [200ms]','FontSize',20)
        xlabel('Intensity Log [mW/mm^2]','FontSize',20)
        legend('Response','Spontaneous Activity');
        title(['Unit ',num2str(i)]);
        ylim([0 max(Data.ProstheticIntensityCount{i})+0.5])
        %   Latency
        if length(find(~isnan(Data.ProstheticLatency{i}))) > 1
        Data.LatencyCurve{i} = figure();
        plot(Data.ProstheticIntensity,Data.ProstheticLatency{i})
        xticks(linspace(0,round(max(Data.ProstheticIntensity),1),5))
        ylabel('Latency [ms]','FontSize',20)
        xlabel('Intensity [mW/mm^2]','FontSize',20)
        legend('Latency of Response');
        title(['Unit ',num2str(i)]);
        ylim([min(Data.ProstheticLatency{i})-10 max(Data.ProstheticLatency{i})+10])
        end
        end
end
end