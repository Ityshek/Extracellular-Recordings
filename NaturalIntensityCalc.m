function [Data] = NaturalIntensityCalc(Data,PSTHbinsize,fname,window,PSTH,sampling_freq,SponFlag2)
%% Natural Intensity Response Curve
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if answer{1} == '1'
    Data.NaturalIntensity = []; Data.NaturalIntensityResponse = cell(1,length(Data.Clusters));   Data.NaturalLatency = cell(1,length(Data.Clusters)); Data.NaturalIntensityCount = cell(1,length(Data.Clusters));
    Data.StimThresh = cell(1,length(Data.Clusters));
end

% Calculation
% Amplitude to Intensity Conversion
% B = [10,20,30,40,50,60,70,80,90,100;0.55,1.2,2.9,5.6,7.9,11.4,18.6,25.6,34.6,42.7]; % 1st row are current in Amp. Second row are intensity in CD/m^2.
% a = [strfind(fname,'Hz')+1,strfind(fname,'amp')];
% Amp = fname(a(1)+2:a(2)-1);
%IntensityTable = [0.55,1.2,2.9,5.6,7.9,11.4,18.6,25.6,34.6,42.7];
if ~isnan(regexp(fname,'_[0123456789.]*CD'))
    [Int1,Int2] = regexp(fname,'_[0123456789.]*CD');
    StimType = 1;
    Data.NaturalIntensity = [Data.NaturalIntensity;str2num(fname(Int1+1:Int2-2))]; % Find Stim Intensity in CD/m^2.
elseif regexp(fname,'_[0123456789.]*nW')~= []
    [Int1,Int2] = regexp(fname,'_\w*Sc');
    StimType = 2;
    Data.NaturalIntensity = [Data.NaturalIntensity;str2num(fname(Int1+1:Int2-2))]; % Find Stim Intensity in micro-watts/mm^2.
end
% Response Calculation
for i=1:length(Data.Clusters)
    ResponseWindow = [(window(1)/10)+2:(window(2)/10)+2]; % Define time window for Prosthetic response (10-210ms post trigger). add 2 bins for -10 and 0 bins in PSTH.
    Data.NaturalIntensityResponse{i} = [Data.NaturalIntensityResponse{i};round(max(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow)),2)];
    Data.NaturalIntensityCount{i} = [Data.NaturalIntensityCount{i};Data.SpikeCount{i}];
end

% Stimulation Threshold
for i=1:length(Data.Clusters)
    %     [RS,RM] = std(Data.Psth_sort{Data.Clusters(i)},'omitnan'); % Calc std and mean of total acticvity
    [Val,Ind]= max(Data.Psth_sort{Data.Clusters(i)}((window(1)/10)+3:(window(2)/10)+2));
    R = sum(Data.Psth_sort{Data.Clusters(i)}(Ind+1:Ind+4));
    count = 1;
    Skips = (window(2)-window(1))/10;
    for k = window(2)/10+1+Skips:Skips:length(Data.Psth_sort{Data.Clusters(i)})
        temp = full(Data.Rast_sort{Data.Clusters(i)}(:,(k-Skips)*0.01*sampling_freq:(k)*0.01*sampling_freq));
        %SponBinned(count) = sum(Data.Psth_sort{Data.Clusters(i)}(k:k+3));
        SponBinned(count) = sum(sum(temp))/size(temp,1);    
        count = count+1;
    end
    if SponFlag2 ~=1
    RM = Data.Spon{1}(end); RS = Data.SponStd{1}(end);
    else
    RM = mean(SponBinned);   RS = std(SponBinned);
    end
    if Data.SpikeCount{Data.Clusters(i)} > RM+2*RS % Look if R window is bigger then 95% of total activity
        Data.StimThresh{i} = [Data.StimThresh{i};1];
        % Calc Latency (Only if responsive)
        Latency = (find(PSTH((window(1)/10)+2:(window(2)/10)+2)...
            == max(PSTH((window(1)/10)+2:(window(2)/10)+2)),1))*PSTHbinsize;
        Data.NaturalLatency{i} = [Data.NaturalLatency{i};Latency];
            else
        Data.StimThresh{i} = [Data.StimThresh{i};0];
        Data.NaturalLatency{i} = [Data.NaturalLatency{i};nan];
    end
end

% Prosthetic Response Latency
% % for i=1:length(Data.Clusters)
% %     V = min(find(round(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow),2) == Data.NaturalIntensityResponse{i}(end))+min(ResponseWindow)-1);
% %     Data.NaturalLatency{i} = [Data.NaturalLatency{i};(V-2)*PSTHbinsize]; % Calculate latency of response in mSec.
% % end

if answer{2} =='1'
    % Plotting
    %   Spiking Rate
    for i=1:length(Data.Clusters)
        Data.IntensityCurve{i} = figure();
        spon = ones(1,length(Data.NaturalIntensityResponse{i}))*(mean(Data.Spon{i}));
%         semilogx([Data.NaturalIntensity],[Data.NaturalIntensityResponse{i}],'-',Data.NaturalIntensity,spon,'--')
        %xticks(linspace(0,round(max(Data.NaturalIntensity),1),5))
        xticks([0;Data.NaturalIntensity])
        ylabel('Spiking Rate[Hz]','FontSize',20)
        if StimType == 1
        xlabel('Intensity log [CD/m^2]','FontSize',20)
        else
            xlabel('Intensity log [nW/mm^2]','FontSize',20)
        end
        legend('Intensity Response')%,'Spontaneous Activity');
        title(['Unit ',num2str(i)]);
        ylim([max(Data.Spon{i})-5 max(Data.NaturalIntensityResponse{i})+5])
        xlim([0 max(Data.NaturalIntensity)])
        % Spike Count
        Data.IntensityCurve{i} = figure();
        spon = ones(1,length(Data.NaturalIntensityResponse{i}))*(mean(Data.Spon{i}));
%         semilogx([Data.NaturalIntensity],[Data.NaturalIntensityCount{i}],'-',Data.NaturalIntensity,spon,'--')
        %xticks(linspace(0,round(max(Data.NaturalIntensity),0),6))
        xticks([0;Data.NaturalIntensity])
        ylabel('Spike Count','FontSize',20)
        if StimType == 1
        xlabel('Intensity log [CD/m^2]','FontSize',20)
        else
            xlabel('Intensity log [nW/mm^2]','FontSize',20)
        end
        %legend('Intensity Response')%,'Spontaneous Activity');
        title(['Unit ',num2str(i)]);
        ylim([0 max(Data.NaturalIntensityCount{i})+0.1])
        xlim([0 max(Data.NaturalIntensity)])
        %   Latency
        Data.LatencyCurve{i} = figure();
        plot(Data.NaturalIntensity,Data.NaturalLatency{i})
        xticks(linspace(0,round(max(Data.NaturalIntensity),1),5))
        ylabel('Latency[ms]','FontSize',20)
        if StimType == 1
        xlabel('Intensity log [CD/m^2]','FontSize',20)
        else
            xlabel('Intensity log [nW/mm^2]','FontSize',20)
        end
        legend('Latency of Response');
        title(['Unit ',num2str(i)]);
        ylim([min(Data.NaturalLatency{i})-10 max(Data.NaturalLatency{i})+10])
    end
end
end
