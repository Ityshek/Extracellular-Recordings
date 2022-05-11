close all;
clearvars -except ProstheticIntensityResponse ProstheticIntensity CPDs CPDResponse;
RastFig = figure; RawFig = figure; PsthFig = figure;
[fname,pathname]=uigetfile('*.mat','Choose data file');
ChannelPosition = [1,9,5,6,14,13,2,10,15,7,12,11,3,4,16,8];
std_Factor=-5;
FigPlotNum = 1;
for c = 1:1 % Change according to the number of recorded channels
% Select One of the two rows below according to the stimulation system 
    
   [raw_data, sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag] =load_data_MultiUnit(c,fname,pathname); % data loader for Screen stimulation
   %[raw_data,sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag,stimulus_times, stimulus_indexes] = SUload_data_Micron(c,fname,pathname); % data loader for Micron stimulation

    if ~channelflag
        [stimulus_times,stimulus_indexes]=find_stim(stim_Data,stim_sampling_rate,sampling_freq,Begin_record); % run only for screen stimulus

        [startIndex,endIndex] = regexp(fname,'_\d*st');
        if startIndex ~= 0
            StimDuration = str2double(fname(startIndex+1:endIndex-2));
        else
            StimDuration = 0;
        end
       
        stimulus_indexes = stimulus_indexes-round(StimDuration/1000*stim_sampling_rate);
        t=[0:length(raw_data)-1]/sampling_freq;        
        Data.thresh(c)=std_Factor*nanstd(double(raw_data));
        %% Plot raw data
        stim=nan(length(raw_data),1);
        stim(stimulus_indexes)=90;
        figure(RawFig);
        subplot(FigPlotNum,FigPlotNum,ChannelPosition(c))
        threshold = Data.thresh(c)*ones(1,length(t));
        plot(t,raw_data,'b',t,stim,'^g',t,threshold,'r')

        % figure settings
        xlabel('Time[Sec]','FontSize',20)
        ylabel('Amplitude[\muV]','FontSize',20)
        title(['Channel ',num2str(c)])
        %         legend('Raw Signal','Trigger','Detected Spikes','Threshold')
        ylim([-400 400])
        xlim([0 max(t)])
        %% Build raster
        outlier = str2num(cell2mat(inputdlg('Set Amplitude threshold to filter out outliers')));
        [Rast,Spike,Av_spike,indx_spike,ind_rast,spike_stim,spike_times]=build_rastRef(stimulus_indexes,stimulus_times,raw_data,sampling_freq,Data.thresh(c),outlier);
        Data.SpikeTimes{c} = spike_times;
        Data.IndxSpike{c} =  indx_spike;
        Data.Spike{c} = Spike;
        figure(RastFig);
        subplot(FigPlotNum,FigPlotNum,ChannelPosition(c))
        spy(Rast)
        x=(set(gca, 'XTickLabel',(linspace(-10*10^-3,mean(diff(stimulus_times)),10))));
        xtickc=linspace(-10*10^-3*sampling_freq,round(mean(diff(stimulus_times(2:end))),2)*sampling_freq,10);
        names= round(((xtickc/sampling_freq)-10^-3),2)*10^3;
        
        % figure settings
        set(gca, 'XTick',  xtickc, 'XTickLabel', names)
        axis square
        xlabel('Time[mSec]','FontSize',20)
        ylabel('Stimulus Repetition','FontSize',20)
        title(['Channel ',num2str(c)])
        xlim([0 length(Rast)])
        
        % Sponteneous Activity in Hz, calculated from recording prior to 1st trigger.
        NumSponSpikes = length(find(spike_times<stimulus_times(1)));
        Spon = NumSponSpikes/stimulus_times(1);
        %% Build PSTH
        figure(PsthFig);
        subplot(FigPlotNum,FigPlotNum,ChannelPosition(c))
        [Psth,binsize_sec]=Build_psth3(Rast,sampling_freq);
        Data.PSTH{c} = Psth;
        t_pst=1000*linspace(-10*10^-3,size(Psth,1)*binsize_sec,length(Psth));
        b = bar(t_pst,Psth);
        hold on
        % figure settings
        xticks([-10:100:size(Psth,1)*binsize_sec*10^3])
        xlabel('Time[mSec]','FontSize',20)
        ylabel('Spiking Rate[Hz]','FontSize',20)
        title(['Channel ',num2str(c)])
        ylim([0 200]);

    end
end
        figure(RawFig);
        x = nan(1,length(t));
        x(indx_spike)=raw_data(indx_spike);
        hold on
        plot(t,x,'*k')
%% Responsive Channels Analysis
%close all;
ActiveChannels = cell2mat(inputdlg('Select Channels for Further Analysis'));
ActiveChannels = str2num(ActiveChannels);
FigPlotNum = round(length(ActiveChannels)/2);
AvgSpkFig = figure;  ClusterResultsFig = figure;
%TIHFig = figure;
col=['r','g','b','m','c','y','k'];
for c=1:length(ActiveChannels)
    AvgsortedSpkFig{c} = figure;
    % TIH
    %     Data.ISI{ActiveChannels(c)} = diff(Data.SpikeTimes{ActiveChannels(c)});
    %     figure(TIHFig)
    %     subplot(2,2,c)
    %     edges = [0:0.1:50];
    %     histogram(Data.ISI{ActiveChannels(c)}*10^3,edges);
    %     xlabel('ISI[mSec]')
    %     ylabel('Count')
    
    % Average Spike
    [Data.AlignedSpikes{ActiveChannels(c)},Data.AverageSpike{ActiveChannels(c)},Aligned_idx]=Align_spikes4(Data.Spike{ActiveChannels(c)},sampling_freq,std_Factor, Data.IndxSpike{ActiveChannels(c)});
    t_spike=((0:length(Data.AverageSpike{ActiveChannels(c)})-1)/sampling_freq)*10^3;
    figure(AvgSpkFig);
    subplot(FigPlotNum,FigPlotNum,c)
    plot(t_spike,Data.AlignedSpikes{ActiveChannels(c)})
    hold on
    plot(t_spike,Data.AverageSpike{ActiveChannels(c)},'k','linewidth',2)
    xlabel('Time[mSec]','FontSize',20)
    ylabel('Amplitude[\muV]','FontSize',20)
    ylim([-100 100]);
    title(['Spike Waveforms Channel ',num2str(ActiveChannels(c))]);
    % PCA + Clustering
    %ClustEvalDB = evalclusters(Data.AlignedSpikes{ActiveChannels(c)},'kmeans','DaviesBouldin','KList',[1:5]);
    %ClustEvalSILL = evalclusters(Data.AlignedSpikes{ActiveChannels(c)},'kmeans','Silhouette','KList',[1:5]);
    %Data.dim{ActiveChannels(c)}=round(mean([ClustEvalSILL.OptimalK ClustEvalDB.OptimalK]));
    Data.dim{ActiveChannels(c)}= 2;
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
    for i=1:Data.dim{ActiveChannels(c)}
        t_sort=((0:size(Data.SortedSpikes{1,ActiveChannels(c)}{i},1)-1)/sampling_freq)*10^3;
        Avg_Sorted_Spikes(:,i) = mean(Data.SortedSpikes{ActiveChannels(c)}{i},2);
        figure(AvgsortedSpkFig{c});
        subplot(1,Data.dim{ActiveChannels(c)},i)
        for u = 1:size(Data.SortedSpikes{ActiveChannels(c)}{i},2)
        plot(t_sort,Data.SortedSpikes{ActiveChannels(c)}{i}(:,u))
        hold on
        end
        plot(t_sort,Avg_Sorted_Spikes(:,i),'k','linewidth',2)
        title(['Sorted waveforms - Cluster ',num2str(i)]);
        hold on
        xlabel('Time[mSec]','FontSize',20)
        ylabel('Amplitude[\muV]','FontSize',20)
        ylim([-200 100]);
    end

    

    %  raw data with sorted spikes
%         stim=nan(length(raw_data),1);
%         stim(stimulus_indexes)=90;
        figure();
        subplot(FigPlotNum,FigPlotNum,ChannelPosition(c))
%       threshold = Data.thresh(c)*ones(1,length(t));
        plot(t,raw_data,'b',t,stim,'^g',t,threshold,'r')

        % figure settings
        xlabel('Time[Sec]','FontSize',20)
        ylabel('Amplitude[\muV]','FontSize',20)
        title(['Channel ',num2str(c)])        
        ylim([-400 400])
        xlim([0 max(t)])  
        hold on

        for i=1:Data.dim{ActiveChannels(c)}  
            ClusteredIdx{i} = nan(1,length(Aligned_idx));
            ClusteredVal{i} = nan(1,length(Aligned_idx));
            for k = 1:length(Aligned_idx)
                if  Data.ClusterIdx{c}(k) == i
                ClusteredIdx{i}(k) = Aligned_idx(k);
                ClusteredVal{i}(k) = raw_data(Aligned_idx(k));
            end
        end
             y{i} = nan(1,length(t));
             y{i}(rmmissing(ClusteredIdx{i})) = rmmissing(ClusteredVal{i});
             plot(t,y{i},['*';col(i)]);
             hold on
        end
 legend('Raw Signal','Trigger','Threshold','Detected Spikes Cluster 1','Detected Spikes Cluster 2')
end
%% Sorted Plots
for c = 1:1
[Rast_sort]=build_rast_sort2(Data.ClusterIdx{ActiveChannels(c)},stimulus_indexes,sampling_freq,ind_rast,Data.dim{ActiveChannels(c)},stimulus_times);
for i=1:Data.dim{ActiveChannels(c)}
    sorted_rast(i)=figure;
    spy(Rast_sort{i})
%     x=(set(gca, 'XTickLabel',(linspace(-1*10^-3,1,10))));
%     xtickc=linspace(1*10^-3*sampling_freq,sampling_freq,5);
%     names= round(((xtickc/sampling_freq)-10^-3),2);
    %    xtick=[1*10^-3:(0.5*10^4)/sampling_freq:1];
    %  xticlables=[0:0.25*sampling_freq:sampling_freq];
    set(gca, 'XTick',  xtickc, 'XTickLabel', names)
    
    axis square
    xlabel('Time[mSec]','FontSize',20)
    ylabel('Stimulus Repetition','FontSize',20)
    title(['Sorted raster ', num2str(i)])

end




%% PSTH Sorted
for i=1:Data.dim{ActiveChannels(c)}
    
    Psth_sort{i}=Build_psth(Rast_sort{i},sampling_freq);
    % Psth2=Build_psth(Rast_sort{2},sampling_freq);
    % Psth3=Build_psth(Rast_sort{3},sampling_freq);
    
    sorted_PSTH(i)=figure;
    bar(t_pst,Psth_sort{i})
    xlabel('Time[mSec]','FontSize',20)
    ylabel('Spiking Rate[Hz]','FontSize',20)
    title(['Sorted PSTH ', num2str(i)])
    ylim([0 100])
    %     file=[fname(1:end-4),'_sortedPSTH','_G',num2str(i)];
%     hgsave(sorted_PSTH(i),[path file(1:end),'.fig' ],'-v7.3');
%     saveas(sorted_PSTH(i),[path file(1:end),'.tif' ]);
end

end
%% Prosthetic Intensity Response Curve
    ProstheticIntensity = []; ProstheticIntensityResponse = [];    
    %% Calculation
        % Amplitude to Intensity Conversion
B = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7;1.1,2.9,4.6,6.501,7.9,9.1,10.1,11.7,13,14.4,15.7,17.3,18.9]; % 1st row are current in Amp. Second row are intensity in mW.
a = [strfind(fname,'Hz')+1,strfind(fname,'amp')];
Amp = fname(a(1)+2:a(2)-1);
ProstheticIntensity = [ProstheticIntensity;B(find(B == str2num(Amp))+1)*2.5/9]; %Convert to Intensity/mm^2 for 40% duty cycle.
        % Response Calculation    
ResponseWindow = 0.06/binsize_sec + 2; % Define time window for Prosthetic response. add 2 bins for -10 and 0 bins in PSTH.
ProstheticIntensityResponse = [ProstheticIntensityResponse;max(Data.PSTH{1}(4:ResponseWindow))]; 
    %% Plotting
figure();
plot(ProstheticIntensity,ProstheticIntensityResponse,'-')
xticks(round(ProstheticIntensity,1))
ylabel('Spiking Rate[Hz]','FontSize',20)
xlabel('Intensity[mW/mm^2]','FontSize',20)

    
    
    
    
%% CPD Selectivity over different data files
CPDs = []; CPDResponse = [];
    %% Calculation
a = [max(strfind(fname,'0_')),strfind(fname,'CPD')];
CPD = fname(a(1):a(2)-1);
CPD(strfind(CPD,'_')) = '.';
CPDs = [CPDs; str2num(CPD)];
ResponseWindow = [1.1/binsize_sec:1.1/binsize_sec+4];
CPDResponse = [CPDResponse; max(max(Data.PSTH{1}(ResponseWindow)))];
    %% Plotting
figure();
plot(CPDs,CPDResponse,'-');
xticks(round(linspace(min(CPDs),max(CPDs),10),2));
ylabel('Spiking Rate[Hz]','FontSize',20);
xlabel('CPD','FontSize',20);
ylim([0 150]);
