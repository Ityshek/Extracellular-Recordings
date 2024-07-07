%% Load & Plot Basic Data
close all;
clearvars -except Data;
%[fname,pathname]=uigetfile('*.mat','Choose data file');
std_Factor=-5;
prompt = { 'Number of Recorded Channels:','Number of the First Recoeded Channel:', 'N rows:','N Column:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1','27','1','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput); RawFig = figure;
FigPlotNum = str2double(answer{1}); count = str2double(answer{2}); nRows = str2num(answer{3}); nColumn = str2num(answer{4});
definptOut = {'150'};
raw_data = cell(str2num(answer{1}));
Dir=uigetdir('*.mat','Select a Folder to Load Raw Data From');
for c = 1:FigPlotNum
    RastFig = figure;  PsthFig = figure;
    chan = {num2str(str2double(answer{2})-1+c)};
    RC = [str2double(cell2mat(inputdlg('Insert the Number of the Recorded Channel',dlgtitle,dims,chan)))];
    [raw,sampling_freq,stim_Data,stim_sampling_rate,Begin_record,stimulus_times,stimulus_indexes,StimIdx,StimFreqs,FreqsOrder,Check] = load_data_ConcateCFFNIR(RC,Dir); % data loader for Micron stimulation
    raw_data{c} = raw;
%   if ~channelflag
        t=[0:length(raw_data{c})-1]/sampling_freq;
        Data.thresh(c)=std_Factor*nanstd(double(raw_data{c}));

        % Plot raw data%
        stim=nan(length(raw_data{c}),1);
        stim(stimulus_indexes(2:end))=150;
        figure(RawFig);
        subplot(nRows,nColumn,c)
        threshold{c} = Data.thresh(c)*ones(1,length(t));
        plot(t,raw_data{c},'b',t,stim,'vg',t,threshold{c},'r')

        % figure settings
        %xlabel('Time[Sec]','FontSize',15)
        %ylabel('Amplitude[\muV]','FontSize',15)
        title(['Channel ',num2str(RC)])
        %         legend('Raw Signal','Trigger','Detected Spikes','Threshold')
        ylim([-200 200])
        xlim([0 max(t)])
        xlabel('Time[Sec]','FontSize',20)
        ylabel('Amplitude[\muV]','FontSize',20)

        % Build raster%
        prompt = {'Insert Upper Thershold Value for spike outliers:'};
        dlgtitle = 'Input'; dims = [1 35];
        outlier = str2num(str2mat(inputdlg(prompt,dlgtitle,dims,definptOut)));
        VisFlag{1} = 3;
        for F = 1:length(FreqsOrder)
        
        
        Idxs = find(StimFreqs(1:end-1) == FreqsOrder(F));
        [Rast,Spike,Av_spike,indx_spike,ind_rast,spike_stim,spike_times]=build_rastRef(stimulus_indexes(Idxs),stimulus_times(Idxs),raw_data{c},sampling_freq,Data.thresh(c),outlier,VisFlag,Check);
        Data.SpikeTimes{c} = spike_times;
        Data.IndxSpike{c} =  indx_spike;
        Data.Spike{c} = Spike;
        figure();
%         figure(RastFig);
%         subplot(1,length(FreqsOrder),F)
        spy(Rast)
        x=(set(gca,'XTickLabel',(linspace(-10*10^-3,mean(diff(stimulus_times))-10*10^-3,5))));
        xtickc=linspace(-10*10^-3*sampling_freq,round(median(diff(stimulus_times(2:end))),2)*sampling_freq...
            -10*10^-3*sampling_freq,5);
        names= round(((xtickc/sampling_freq)-10^-3),2)*10^3;
        %xlim([-440 43560]);
        % figure settings
        set(gca, 'XTick',  xtickc, 'XTickLabel', names)
        axis square
        xlabel('Time [ms]','FontSize',10)
        ylabel('Stimulus Repetition','FontSize',10)
        title(['Stim Freq:',num2str(FreqsOrder(F))])




        % Build PSTH%
        CountWindow = [2 22];
        %figure(PsthFig);
        Pfig{F} = figure();
        ax = axes();
        %subplot(1,length(FreqsOrder),F)
        [Psth,binsize_sec,smoothed_Psth,SpikeCount,Label]=Build_psth3(Rast,sampling_freq,CountWindow);
        Data.PSTH{c}{F} = Psth;
        StimTime = round(mean(diff(stimulus_times)),2);
        t_pst=1000*linspace(-10*10^-3,size(Psth,2)*binsize_sec-10*10^-3,length(Psth));
        b = bar(t_pst,Psth);
        hold on
        plot(t_pst,smoothed_Psth,'linewidth',2)
        % figure settings
        %xticklabels(linspace(0,StimTime*1000,5))
        %xticklabels([-10])
        axis square
        xlabel('Time [ms]')
        ylabel(Label)
        title(['Stim Freq:',num2str(FreqsOrder(F)),' | CH:',num2str(RC)])
        ylim([0 1]);  xlim([-10 1900]);
        %       ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
        %       ax.Box = 'off'; ax.Color = "none";
%     end
        end    
count = count+1;
    figure(RawFig);
    x = nan(1,length(t));
    x(indx_spike)=raw_data{c}(indx_spike);
    hold on
    plot(t,x,'*k')

    % Calculate SNR
%     SponIdxSpikes = indx_spike(find(indx_spike<stimulus_indexes(2)));
%     SponSpike = Spike(1:length(SponIdxSpikes));
%     [Data.SNR{c}] = SNRCalc(SponSpike,raw_data{c}(1:stimulus_indexes(2)),Data.thresh(c),SponIdxSpikes);
%% PCA & Clustering
%close all;
definput2 = "1"; prompt = {'Select Channels for Further Analysis'};
ActiveChannels = cell2mat(inputdlg(prompt,dlgtitle,dims,definput2));
ActiveChannels = str2num(ActiveChannels);
FigPlotNum = length(ActiveChannels);
AvgSpkFig = figure;  ClusterResultsFig = figure;
col=['r','g','m','c','y','k'];
for c=1:length(ActiveChannels)
    % Average Spike
    [Data.AlignedSpikes{ActiveChannels(c)},Data.AverageSpike{ActiveChannels(c)},Aligned_idx]=Align_spikes4(Data.Spike{ActiveChannels(c)},sampling_freq,std_Factor, Data.IndxSpike{ActiveChannels(c)});
    t_spike=((0:length(Data.AverageSpike{ActiveChannels(c)})-1)/sampling_freq)*10^3;
    figure(AvgSpkFig);
    subplot(FigPlotNum,FigPlotNum,c)
    plot(t_spike,Data.AlignedSpikes{ActiveChannels(c)})
    hold on
    plot(t_spike,Data.AverageSpike{ActiveChannels(c)},'k','linewidth',2)
    xlabel('Time [ms]','FontSize',20)
    ylabel('Amplitude [\muV]','FontSize',20)
    ylim([-150 150]);
    title(['Spike Waveforms Channel ',num2str(ActiveChannels(c))]);

    % PCA + Clustering
    ClustEvalCH = evalclusters(Data.AlignedSpikes{ActiveChannels(c)},'kmeans','CalinskiHarabasz','KList',[1:5]);
    ClustEvalG = evalclusters(Data.AlignedSpikes{ActiveChannels(c)},'kmeans','gap','KList',[1:5]);
    Data.dim{ActiveChannels(c)}=min([ClustEvalG.OptimalK ClustEvalCH.OptimalK]);
    %Data.dim{ActiveChannels(c)}= 2;

    figure(ClusterResultsFig)
    subplot(FigPlotNum,FigPlotNum,c)
    [Data.ClusterIdx{ActiveChannels(c)},C,score]=PCA_Analysis5(Data.AlignedSpikes{ActiveChannels(c)},Data.dim{ActiveChannels(c)});
    title(['Cluster Results Channel',num2str(ActiveChannels(c))]);
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    hold on
    % Sorted Waveforms
    Data.SortedSpikes{ActiveChannels(c)}=sort_spikes3(Data.ClusterIdx{ActiveChannels(c)},Data.AlignedSpikes{ActiveChannels(c)},Data.dim{ActiveChannels(c)});
    color = {[1.00,0.42,0.42],[0.42,1,0.42],'b'};
    color2 = {'r','g','b'};
    AvgsortedSpkFig{c} = figure;
    for i=1:Data.dim{ActiveChannels(c)}
        t_sort=((0:size(Data.SortedSpikes{1,ActiveChannels(c)}{i},1)-1)/sampling_freq)*10^3;
        Avg_Sorted_Spikes(:,i) = mean(Data.SortedSpikes{ActiveChannels(c)}{i},2);
        figure(AvgsortedSpkFig{c});
        subplot(1,Data.dim{ActiveChannels(c)},i)
        for u = 1:size(Data.SortedSpikes{ActiveChannels(c)}{i},2)
            plot(t_sort,Data.SortedSpikes{ActiveChannels(c)}{i}(:,u),'Color',color{i})
            hold on
        end
        plot(t_sort,Avg_Sorted_Spikes(:,i),'k','linewidth',10,'Color',color2{i})
        title(['Sorted waveforms - Cluster ',num2str(i)]);
        hold on
        xlabel('Time[mSec]','FontSize',20)
        ylabel('Amplitude[\muV]','FontSize',20)
        ylim([-250 250]); xlim([1 2.5]);
    end



    %  raw data with sorted spikes
    %         stim=nan(length(raw_data),1);
    %         stim(stimulus_indexes)=90;
    figure();
    subplot(FigPlotNum,FigPlotNum,c)
    %       threshold = Data.thresh(c)*ones(1,length(t));
    plot(t,raw_data{c},'b',t,stim,'vg',t,threshold{c},'r')

    % figure settings
    xlabel('Time [sec]','FontSize',20)
    ylabel('Amplitude[\muV]','FontSize',20)
    title(['Channel ',num2str(c+16)])
    ylim([-250 250])
    xlim([0 max(t)])
    hold on

    for i=1:Data.dim{ActiveChannels(c)}
        Data.ClusteredspikeIdx{i} = nan(1,length(Aligned_idx));
        Data.ClusteredVal{i} = nan(1,length(Aligned_idx));
        for k = 1:length(Aligned_idx)
            if  Data.ClusterIdx{ActiveChannels(c)}(k) == i
                Data.ClusteredspikeIdx{i}(k) = Aligned_idx(k);
                Data.ClusteredVal{i}(k) = raw_data{c}(Aligned_idx(k));
            end
        end
        y{i} = nan(1,length(t));
        y{i}(rmmissing(Data.ClusteredspikeIdx{i})) = rmmissing(Data.ClusteredVal{i});
        plot(t,y{i},['*';col(i)]);
        hold on
    end
    legend('Raw Signal','Trigger','Threshold','Detected Spikes Cluster 1','Detected Spikes Cluster 2')
    %   TIH
    TIHFig = figure;
    for i=1:Data.dim{ActiveChannels(c)}
        Data.ISI{ActiveChannels(c)}{i} = (diff(rmmissing(Data.ClusteredspikeIdx{i}))/sampling_freq)*10^3; % Save the ISIs in ms for every cluster
        figure(TIHFig)
        subplot(2,2,i)
        edges = [0:0.5:100];
        histogram(Data.ISI{ActiveChannels(c)}{i}...
            ,edges,'Normalization','pdf');
        xlabel('ISI [ms]')
        ylabel('Probability')
        ylim([0 0.3])
    end
    %CorrPlot = CorrFunc(t,Data.ClusteredspikeIdx,Data.dim{ActiveChannels(c)},sampling_freq);
end
%% Select Relevant Clusters
Data.Clusters = [str2num(cell2mat(inputdlg('Insert the Number of the Relevant Clusters')))]; % For multi file analysis, insert cluster numbers according of the unit's order from previous files.
% if ~isfield(Data,'StimThresh')
%     Data.Rast_sort = {}; Data.Psth_sort = {}; Data.ind_rast_sort ={};
%     Data.Rast_sort{end+1} = Rast_sort{Data.Clusters};
%     Data.Psth_sort{end+1} = Psth_sort(Data.Clusters,:);
%     Data.ind_rast_sort{end+1} = ind_rast_sort{Data.Clusters};
% else
%     Data.Rast_sort{end+1} = Rast_sort{Data.Clusters};
%     Data.Psth_sort{end+1} = Psth_sort(Data.Clusters,:);
%     Data.ind_rast_sort{end+1} = ind_rast_sort{Data.Clusters};
% end
%% Sponteneous Activity calculated from recording prior to 1st trigger.
% Run next line Only for first trial of each experiment
% if ~isfield(Data,'StimThresh')
%     Data.Spon = {}; Data.SponStd = {};
% end
% % SponFlag = cell2mat(inputdlg('Choose Natural / Prosthetic (1 = N | 0 = P)'));
% % if SponFlag == '0'
% SponWindow = (window(2)-window(1))*0.001;
% LastSponSpike = max(find(spike_times<stimulus_times(1)));
% 
% [Spon,SponStd] = SponCalc(Data,sampling_freq,stimulus_indexes,window);
% Data.Spon{end+1} = Spon{1}; Data.SponStd{end+1} = SponStd{1};
% %
% % for i=1:length(Data.Clusters)
% %     SponCount = length(find(Data.ClusteredspikeIdx{i} < stimulus_indexes(1)));
% %     SponSpikeTrain = zeros(1,int32(spike_times(LastSponSpike)*sampling_freq));
% %     SponSpikeTrain(int32(spike_times(find(Data.ClusteredspikeIdx{i}...
% %         <=spike_times(LastSponSpike)*sampling_freq))*sampling_freq)) = 1; % Find all spikes of the current cluster before stimulus onset
% %
% %     NumBins = round(spike_times(LastSponSpike)/SponWindow);
% %     BinSizeSpon = int32(size(SponSpikeTrain,2)/NumBins); countSpon = 1;
% %     for k=BinSizeSpon:BinSizeSpon:size(SponSpikeTrain,2) - BinSizeSpon
% %         SponBinned(countSpon) = sum(SponSpikeTrain(k-BinSizeSpon+1:k+BinSizeSpon));  % Calc num of Spon spikes in 200ms window (should be same window as Spike Count Var).
% %         countSpon = countSpon+1;
% %     end
% %     Data.Spon{i} = [Data.Spon{i};mean(SponBinned)];
% %     Data.SponStd{i} = [Data.SponStd{i};std(SponBinned)];
% %
% % end
% % else
% %     for i=1:length(Data.Clusters)
% %         for k = 1:4:length(Data.Psth_sort{Data.Clusters(i)}(window(2)/10:end))
% %             SponBinned(k) = sum(Data.Psth_sort{Data.Clusters(i)}(k:k+4));
% %             Data.Spon{i} = [Data.Spon{i};mean(SponBinned)];
% %             Data.SponStd{i} = [Data.SponStd{i};std(SponBinned)];
% %         end
% %     end
% % end
%% Calc Spon by Spikes/sec
Spon1 = length(find(Data.ClusteredspikeIdx{Data.Clusters}<stimulus_indexes(1)/4));
Spon2 = length(find(Data.ClusteredspikeIdx{Data.Clusters}<stimulus_indexes(1)/4*2))-Spon1;
Spon3 = length(find(Data.ClusteredspikeIdx{Data.Clusters}<stimulus_indexes(1)/4*3))-(Spon1+Spon2);
Spon4 = length(find(Data.ClusteredspikeIdx{Data.Clusters}<stimulus_indexes(1)))-(Spon1+Spon2+Spon3);
SponTime = stimulus_times(1)/4;
Spon1 = Spon1/SponTime; Spon2 = Spon2/SponTime; Spon3 = Spon3/SponTime; Spon4 = Spon4/SponTime;
SponStd = std([Spon1 Spon2 Spon3 Spon4]); SponMean = mean([Spon1 Spon2 Spon3 Spon4]);
% Convert to count by Psth binsize
SponFactor = 1/(binsize_sec*5);
SponNew = (SponMean/SponFactor)+(2*(SponStd/SponFactor));
% Draw Baseline in PSTH figures
Baseline = ones(1,length(smoothed_Psth))*SponNew;
for i=1:length(Pfig)
figure(Pfig{i});
hold on
plot(t_pst,Baseline,'Color','k','LineWidth',2)
end
end
% figure(RawFig);
% legend('Raw Signal','Trigger',['Threshold (SNR = ',num2str(Data.SNR{c}),')'],'Detected Spikes');
%% Save Results
%save('C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Results\Data\CFF\CFF_2023.06.12_CH_22_27_Loc4.mat')