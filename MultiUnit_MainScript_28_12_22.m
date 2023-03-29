%% Load & Plot Basic Data
close all;
clearvars -except Data;
RastFig = figure; RawFig = figure; PsthFig = figure;
[fname,pathname]=uigetfile('*.mat','Choose data file');
ChannelPosition = [1,9,5,6,14,13,2,10,15,7,12,11,3,4,16,8];
std_Factor=-4.5;
prompt = {'Select Projection System (Screen = 1, Micron =2)', 'Number of Recorded Channels:','Number of the First Recoeded Channel:', 'N rows:','N Column:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1','1','1','1','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
FigPlotNum = str2double(answer{2}); count = str2double(answer{3}); nRows = str2num(answer{4}); nColumn = str2num(answer{5});
if answer{1} == '2'
prompt = {'Insert stimulus type (1 for NIR/2 for Vis)','Is the trigger signal distorted? (1=YES, 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput1 = {'1','0'};  
answer1 = inputdlg(prompt,dlgtitle,dims,definput1);
end
for c = 1:FigPlotNum
    chan = {num2str(str2double(definput{3})-1+c)};
    RC = [str2double(cell2mat(inputdlg('Insert the Number of the Recorded Channel',dlgtitle,dims,chan)))];

    % SELECT ONE OF THE TWO ROWS BELOW ACCORDING TO THE STIMULATION SYSTEM
    if answer{1} == '1'
    [raw_data{c}, sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag] =load_data_MultiUnit(RC,fname,pathname); % data loader for Screen stimulation
    for i = 1:length(raw_data{1}); if abs(raw_data{1}(i)) > 500; raw_data{1}(i) = mean(raw_data{1}); end; end
    [stimulus_times,stimulus_indexes]=find_stim(stim_Data,stim_sampling_rate,sampling_freq,Begin_record); % run only for screen stimulus
    VisFlag = cell(1);
    else
    [raw_data{c},sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag,stimulus_times, stimulus_indexes,VisFlag] = SUload_data_Micron(RC,fname,pathname,answer1); % data loader for Micron stimulation
    end

    if ~channelflag


        [startIndex,endIndex] = regexp(fname,'_\d*ms');
        if startIndex ~= 0
            StimDuration = str2double(fname(startIndex+1:endIndex-2));
        else
            StimDuration = 0; 
        end
        %stimulus_indexes = stimulus_indexes+round(StimDuration/1000*stim_sampling_rate);
        %stimulus_times = stimulus_times - (1-StimDuration*10^-3); 
        t=[0:length(raw_data{c})-1]/sampling_freq;
        Data.thresh(c)=std_Factor*nanstd(double(raw_data{c}));

        % Plot raw data%
        stim=nan(length(raw_data{c}),1);
        stim(stimulus_indexes(2:end))=90;
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

        % Build raster%
        prompt = {'Insert Upper Thershold Value for spike outliers:'};
        definpt = {'200'}; dlgtitle = 'Input'; dims = [1 35];
        outlier = str2num(str2mat(inputdlg(prompt,dlgtitle,dims,definpt)));
        [Rast,Spike,Av_spike,indx_spike,ind_rast,spike_stim,spike_times]=build_rastRef(stimulus_indexes,stimulus_times,raw_data{c},sampling_freq,Data.thresh(c),outlier,VisFlag);
        Data.SpikeTimes{c} = spike_times;
        Data.IndxSpike{c} =  indx_spike;
        Data.Spike{c} = Spike;
        figure(RastFig);
        subplot(nRows,nColumn,c)
        spy(Rast)
        x=(set(gca, 'XTickLabel',(linspace(-10*10^-3,mean(diff(stimulus_times))-10*10^-3,5))));
        xtickc=linspace(-10*10^-3*sampling_freq,round(mean(diff(stimulus_times(2:end))),2)*sampling_freq...
            -10*10^-3*sampling_freq,6);
        names= round(((xtickc/sampling_freq)-10^-3),2)*10^3;
        xlim([-440 43560]);
        % figure settings
        set(gca, 'XTick',  xtickc, 'XTickLabel', names)
        axis square
        xlabel('Time [ms]','FontSize',15)
        ylabel('Stimulus Repetition','FontSize',15)
        title(['Channel ',num2str(RC)])
        xlim([0 length(Rast)])



        % Build PSTH%
        CountWindow = [2 22];
        figure(PsthFig);
        subplot(nRows,nColumn,c)
        [Psth,binsize_sec,smoothed_Psth,SpikeCount,Label]=Build_psth3(Rast,sampling_freq,CountWindow);
        Data.PSTH{c} = Psth;
        StimTime = round(mean(diff(stimulus_times)),2);
        t_pst=1000*linspace(-10*10^-3,size(Psth,2)*binsize_sec-10*10^-3,length(Psth));
        b = bar(t_pst,Psth);
        hold on
        plot(t_pst,smoothed_Psth,'linewidth',2)
        % figure settings
        %xticklabels(linspace(0,StimTime*1000,5))
        %xticklabels([-10])
        xlabel('Time [ms]','FontSize',15)
        ylabel(Label,'FontSize',15)
        title(['Channel ',num2str(RC)])
        %ylim([0 150]);

    end
    count = count+1;
    figure(RawFig);
    x = nan(1,length(t));
    x(indx_spike)=raw_data{c}(indx_spike);
    hold on
    plot(t,x,'*k')
    
% Calculate SNR
    [Data.SNR{c}] = SNRCalc(Spike,raw_data{c}(1:stimulus_indexes(2)),Data.thresh(c));
    
end
figure(RawFig);
legend('Raw Signal','Trigger',['Threshold (SNR = ',num2str(Data.SNR{c}),')'],'Detected Spikes');


%% PCA & Clustering
%close all;
ActiveChannels = cell2mat(inputdlg('Select Channels for Further Analysis'));
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
    xlabel('Time [mSec]','FontSize',20)
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
        ylim([-250 250]); xlim([1 2.5]); xticklabels(['0';'0.5';'1';'1.5'])
    end



    %  raw data with sorted spikes
    %         stim=nan(length(raw_data),1);
    %         stim(stimulus_indexes)=90;
    figure();
    subplot(FigPlotNum,FigPlotNum,c)
    %       threshold = Data.thresh(c)*ones(1,length(t));
    plot(t,raw_data{c},'b',t,stim,'vg',t,threshold{c},'r')

    % figure settings
    xlabel('Time[Sec]','FontSize',20)
    ylabel('Amplitude[\muV]','FontSize',20)
    title(['Channel ',num2str(c+16)])
    ylim([-400 400])
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

%% Sorted Plots
for c = 1:length(ActiveChannels)
    [Data.Rast_sort]=build_rast_sort2(Data.ClusterIdx{ActiveChannels(c)},stimulus_indexes,sampling_freq,ind_rast,Data.dim{ActiveChannels(c)},stimulus_times);
    for i=1:Data.dim{ActiveChannels(c)}
        sorted_rast(i)=figure;
        spy(Data.Rast_sort{i})
        %     x=(set(gca, 'XTickLabel',(linspace(-1*10^-3,1,10))));
        %     xtickc=linspace(1*10^-3*sampling_freq,sampling_freq,5);
        %     names= round(((xtickc/sampling_freq)-10^-3),2);
        %    xtick=[1*10^-3:(0.5*10^4)/sampling_freq:1];
        %  xticlables=[0:0.25*sampling_freq:sampling_freq];
        set(gca, 'XTick',  xtickc, 'XTickLabel', names)

        axis square
        xlabel('Time [ms]','FontSize',20)
        ylabel('Stimulus Repetition','FontSize',20)
        title(['Sorted raster ', num2str(i)])

    end




    % PSTH Sorted
    for i=1:Data.dim{ActiveChannels(c)}
        
        [PSTH,binsize_sec,Data.Psth_sort{i},Data.SpikeCount{i},label,window]=Build_psth(Data.Rast_sort{i},sampling_freq);
        t_pst=1000*linspace(-10*10^-3,size(Data.Psth_sort{i},2)*binsize_sec-10*10^-3,...
            length(Data.Psth_sort{i}));
        sorted_PSTH(i)=figure;
        bar(t_pst,PSTH)
        hold on
        plot(t_pst,Data.Psth_sort{i},'linewidth',2)
        xlabel('Time [ms]','FontSize',20)
        ylabel(Label,'FontSize',20)
        title(['Sorted PSTH ', num2str(i)])
        %ylim([0 0.5])
        %     file=[fname(1:end-4),'_sortedPSTH','_G',num2str(i)];
        %     hgsave(sorted_PSTH(i),[path file(1:end),'.fig' ],'-v7.3');
        %     saveas(sorted_PSTH(i),[path file(1:end),'.tif' ]);
    end

end

%% Select Relevant Clusters
Data.Clusters = [str2num(cell2mat(inputdlg('Insert the Number of the Relevant Clusters')))]; % For multi file analysis, insert cluster numbers according of the unit's order from previous files.
%% Sponteneous Activity calculated from recording prior to 1st trigger.
% Run next line Only for first trial of each experiment
if ~isfield(Data,'StimThresh')
Data.Spon = cell(1,length(Data.Clusters)); Data.SponStd = cell(1,length(Data.Clusters));
end
% SponFlag = cell2mat(inputdlg('Choose Natural / Prosthetic (1 = N | 0 = P)'));
% if SponFlag == '0'
SponWindow = (window(2)-window(1))*0.001;
LastSponSpike = max(find(spike_times<stimulus_times(1)));
for i=1:length(Data.Clusters)
    SponCount = length(find(Data.ClusteredspikeIdx{i} < stimulus_indexes(1)));
    SponSpikeTrain = zeros(1,int32(spike_times(LastSponSpike)*sampling_freq));
    SponSpikeTrain(int32(spike_times(find(Data.ClusteredspikeIdx{i}...
        <=spike_times(LastSponSpike)*sampling_freq))*sampling_freq)) = 1; % Find all spikes of the current cluster before stimulus onset
   
    NumBins = round(spike_times(LastSponSpike)/SponWindow);
    BinSizeSpon = int32(size(SponSpikeTrain,2)/NumBins); countSpon = 1;
    for k=BinSizeSpon:BinSizeSpon:size(SponSpikeTrain,2) - BinSizeSpon
        SponBinned(countSpon) = sum(SponSpikeTrain(k-BinSizeSpon+1:k+BinSizeSpon));  % Calc num of Spon spikes in 200ms window (should be same window as Spike Count Var).
        countSpon = countSpon+1; 
    end
    Data.Spon{i} = [Data.Spon{i};mean(SponBinned)];
    Data.SponStd{i} = [Data.SponStd{i};std(SponBinned)];
end
% else
%     for i=1:length(Data.Clusters)
%         for k = 1:4:length(Data.Psth_sort{Data.Clusters(i)}(window(2)/10:end))
%             SponBinned(k) = sum(Data.Psth_sort{Data.Clusters(i)}(k:k+4));  
%             Data.Spon{i} = [Data.Spon{i};mean(SponBinned)];
%             Data.SponStd{i} = [Data.SponStd{i};std(SponBinned)];
%         end
%     end
% end

%% Natural Vis Intensity Curve
RateCalcWindow = 10; % PSTH Bin Size in ms 
if StimTime > 0.5
    SponFlag2 = 1;
else 
    SponFlag2 = 0;
end
[Data] = NaturalIntensityCalc(Data,RateCalcWindow,fname,window,PSTH,sampling_freq,SponFlag2);
%% Prosthetic Intensity Response Curve
RateCalcWindow = 10; % PSTH Bin Size in ms
[Data] = ProstheticIntensityCalc(Data,RateCalcWindow,fname,window,PSTH);
%% CPD Selectivity over different data files
RateCalcWindow = 10; % PSTH Bin Size in ms
[Data] = CPDCalc(Data,RateCalcWindow,fname,sampling_freq);
%% Save Data variable
% Change "SavePath" to the desiered folder
SavePath ='C:\Users\Itay\Desktop\Yossi Mandel Lab\Thesis\Data Files\Natural Intensity\SpikeCount\';
stimfreq = num2str(round(1/mean(diff(stimulus_times)))); 
[a,b] = regexp(fname,'_\d*ms'); Data.StimDur = fname(a+1:b-2); 
[a,b] = regexp(pathname,'Data\\\S{10,10}'); date = pathname(a+5:b);
save([SavePath,'NaturalFullFlash_',Data.StimDur,'ms',stimfreq,'Hz',date,'.mat'],"Data");

 %% Write Results into Table & Save Data in a file
% if Data.Clusters > 1
% sheet = table(Data.ProstheticIntensity,Data.ProstheticIntensityResponse{1},Data.ProstheticLatency{1},Data.Spon{1}...
%     ,Data.ProstheticIntensityResponse{2},Data.ProstheticLatency{1},Data.Spon{2},...
%     'VariableNames',{'Intensity[mW/mm^2]','Response[Spikes/Sec] Unit 1','Latency Unit[ms] 1','Baseline Activity[Spikes/Sec] Unit 1'...
%     ,'Response[Spikes/Sec] Unit 2','Latency Unit[ms] 2','Baseline Activity[Spikes/Sec] Unit 2'});
% else
% sheet = table(Data.ProstheticIntensity,Data.ProstheticIntensityResponse{1},Data.ProstheticLatency{1},Data.Spon{1}...
%     ,'VariableNames',{'Intensity[mW/mm^2]','Response[Spikes/Sec] Unit 1','Latency Unit[ms] 1','Baseline Activity[Spikes/Sec] Unit 1'});    
% end
% Path = 'D:\Yossi Mandel Lab\Thesis\Result Tables\';
% Date = pathname(38:47);
% prompt = {'Stimulus Type:','Stim Duration:','Stim Frequency:'};
% definput = {'ProstheticFullFlash','10ms','2Hz'};
% dlgtitle = 'Input';
% dims = [1 35];
% AA = inputdlg(prompt,dlgtitle,dims,definput);
% TableName = [Path,AA{1},'_',AA{2},AA{3},Date,'.xlsx'];
% writetable(sheet,TableName,'AutoFitWidth',1);
% %                                                           %
% Path2 ='D:\Yossi Mandel Lab\Thesis\Data Files\';
% filename = [Path2,AA{1},AA{2},AA{3},Date,'.mat'];
% save(filename,'Data');
% for i=1:length(Data.Clusters)
% savefig(Data.IntensityCurve{i},['D:\Yossi Mandel Lab\SU Data\Raw Data\',Date,'\Results\',AA{1},AA{2},AA{3},'Unit',num2str(i),' Curve']);
% end