close all;
clear all;
std_Factor = -5.5;
NumOrientations = 12;
NumReps = 5;
StimsPerCPD = NumOrientations*NumReps;
col=['r','g','b','m','c','y','k'];
Orientations = [0:30:330];
RC = [str2num(cell2mat(inputdlg('Insert the Numbers of the Recorded Channels')))];
%     PeakCPDFig = figure();
%     RastFig{k} = figure(); OrTuningCurveFig{k} = figure();
%     [startIndex,endIndex] = regexp(fname{k},'_\d*CPD');
%     CPD(k) = str2double(fname{k}(startIndex+1:endIndex-3))/10;
[raw_data, sampling_freq,stim_Data,stim_sampling_rate,Begin_record,stimulus_times,stimulus_indexes,CPDs] =load_data_ConcateMultiUnit(RC);
CPDs = unique(CPDs);
OrientationStimNums{1} = [1:12:StimsPerCPD*length(CPDs)];
OrientationStimNums{1} = reshape(OrientationStimNums{1},NumReps,length(CPDs));
for a = 2:NumOrientations
    OrientationStimNums{a} = OrientationStimNums{1}+a-1;
    OrientationStimNums{a} = reshape(OrientationStimNums{a},NumReps,length(CPDs));
end
% Slot = str2num(cell2mat(inputdlg('Select Slot Number of headstage')));
% if Slot == 2
%     RC = RC-16;
% end
RawDataFig = figure();
IdxISI = StimsPerCPD; % Select the last index of the stimulus times to be included in ISI calculation(because of concatination reasons).
count = 1;
for c = 1:length(RC)
    Data.thresh(RC(c))=std_Factor*nanstd(double(raw_data{RC(c)}));
    prompt = {'Insert Upper Thershold Value for spike outliers:'}; 
    definpt = {'200'}; dlgtitle = 'Input'; dims = [1 35];
    outlier = str2num(str2mat(inputdlg(prompt,dlgtitle,dims,definpt))); 
    % Draw Raw Data Plots
    figure(RawDataFig)
    t=[0:length(raw_data{RC(c)})-1]/sampling_freq;
    stim=nan(length(raw_data{RC(c)}),1);
    stim(stimulus_indexes{RC(c)})=200;
    subplot(ceil(length(RC)/2),ceil(length(RC)/2),count)
    threshold = Data.thresh(RC(c))*ones(1,length(t));
    plot(t,raw_data{RC(c)},'b',t,threshold,'r')
    hold on
     % figure settings
    xlabel('Seconds')
    ylabel('Amplitude[\muV]')
    title(['Channel ',num2str(RC(c))])
    ylim([-450 450])
    xlim([0 max(t)])
    [Rast,Data.Spike{RC(c)},Av_spike,Data.IndxSpike{RC(c)},ind_rast,spike_stim,spike_times]=build_rastRef28_9(stimulus_indexes{RC(c)},raw_data{RC(c)},sampling_freq,Data.thresh(RC(c)),IdxISI,outlier);    
    x = nan(1,length(t));
    x(Data.IndxSpike{RC(c)})=raw_data{RC(c)}(Data.IndxSpike{RC(c)});
    plot(t,x,'*k','MarkerSize',1) % Plot the event times in black
    hold on
    plot(t,stim,'xg','MarkerSize',3) % Plot the trigger times in green
   count = count + 1;
end
%% Spike Sorting
AvgSpkFig = figure(); ClusterResultsFig = figure();
BlankScreenSec = 1;
prompt = {'Select Channels for Further Analysis:'};
AC = str2num(str2mat(inputdlg(prompt)));
prompt = {'Select Row and Column counts for plots:'};
RowsColumns = str2num(str2mat(inputdlg(prompt)));
Count = 1;
for c=1:length(AC)
    CountUnit = 1;
    % Average Spike Waveform
    AvgsortedSpkFig{c} = figure();
    [Data.AlignedSpikes{AC(c)},Data.AverageSpike{AC(c)},Aligned_idx]=Align_spikes4(Data.Spike{AC(c)},sampling_freq,std_Factor, Data.IndxSpike{AC(c)});
    t_spike=((0:length(Data.AverageSpike{AC(c)})-1)/sampling_freq)*10^3;
    figure(AvgSpkFig);
    subplot(RowsColumns(1),RowsColumns(2),Count)
    plot(t_spike,Data.AlignedSpikes{AC(c)})
    hold on
    plot(t_spike,Data.AverageSpike{AC(c)},'k','linewidth',2)
    xlabel('Time[mSec]')
    ylabel('Amplitude[\muV]')
    ylim([-450 450]);
    title(['Spike Waveforms - Channel ',num2str(AC(c))]);


    
    % PCA + Clustering
    ClustEvalDB = evalclusters(Data.AlignedSpikes{AC(c)},'kmeans','CalinskiHarabasz','KList',[1:5]);
    ClustEvalSILL = evalclusters(Data.AlignedSpikes{AC(c)},'kmeans','gap','KList',[1:5]);
    Data.dim{AC(c)}=floor(mean([ClustEvalSILL.OptimalK ClustEvalDB.OptimalK]));
    %Data.dim{RC(c)} =1;
    figure(ClusterResultsFig)
    subplot(RowsColumns(1),RowsColumns(2),Count)
    %[Data.ClusterIdx{AC(c)},C,score,Data.AlignedSpikes{AC(c)}]=PCA_Analysis6(Data.AlignedSpikes{AC(c)},Data.dim{AC(c)},outlier);
    [Data.ClusterIdx{AC(c)},C,score]=PCA_Analysis5(Data.AlignedSpikes{AC(c)},Data.dim{AC(c)});

    title(['Cluster Results - Channel ',num2str(AC(c))]);


        
    % Sorted Waveforms
    Data.SortedSpikes{AC(c)}=sort_spikes3(Data.ClusterIdx{AC(c)},Data.AlignedSpikes{AC(c)},Data.dim{AC(c)});
    for i=1:Data.dim{AC(c)}
        t_sort=((0:size(Data.SortedSpikes{1,AC(c)}{i},1)-1)/sampling_freq)*10^3;
        Avg_Sorted_Spikes(:,i) = mean(Data.SortedSpikes{AC(c)}{i},2);
        figure(AvgsortedSpkFig{c});
        subplot(2,Data.dim{AC(c)},CountUnit)
        plot(t_sort,Data.SortedSpikes{AC(c)}{i})
        hold on
        plot(t_sort,Avg_Sorted_Spikes(:,i),'k','linewidth',3)
        title(['sorted waveforms - Channel: ',num2str(AC(c)),' Unit: ',num2str(i)]);
        xlabel('Time[mSec]')
        ylabel('Amplitude[\muV]')
        BottomLim = min(min(Data.SortedSpikes{AC(c)}{i})) - 10;
        TopLim = max(max(Data.SortedSpikes{AC(c)}{i})) + 10;
        ylim([BottomLim TopLim]);

        CountUnit = CountUnit+1;
    end

    %% Rasters & Orientation Tuning Curve
    for d=1:Data.dim{AC(c)}
        RastOrCPDFig{c}{d} = figure('Name',['Channel: ',num2str(AC(c)),'Unit: ',num2str(d)]);
        OrResponsePerCPDFig{c}{d} = figure('Name',['Channel: ',num2str(AC(c)),'Unit: ',num2str(d)]); %ResponsePerCPDFig{c}{d} = figure();
        count = 1;
        for p=1:length(CPDs)
            %[Rast,Data.Spike{AC(c)},Av_spike,Data.IndxSpike{AC(c)},ind_rast,spike_stim,spike_times]=build_rastRef28_9(stimulus_indexes{AC(c)}((p-1)*StimsPerCPD+1:p*StimsPerCPD),raw_data{AC(c)},sampling_freq,Data.thresh(c),IdxISI);
            [Data.SortedRasters{AC(c)}]=build_rast_sort4(Data.ClusterIdx{AC(c)},stimulus_indexes{AC(c)}((p-1)*StimsPerCPD+1:p*StimsPerCPD),sampling_freq,Data.dim{AC(c)},Aligned_idx,stimulus_times{AC(c)}((p-1)*StimsPerCPD+1:p*StimsPerCPD));
            %StimTimeForPSTHCalc = mean(diff(stimulus_times{AC(c)}((p-1)*StimsPerCPD+1:p*StimsPerCPD)))-2;
            for i=1:NumOrientations
                figure(RastOrCPDFig{c}{d})
                subplot(3,4,count)
%                subplot(length(CPDs),NumOrientations,count)
                spy(Data.SortedRasters{AC(c)}{d}(OrientationStimNums{i}(:,1),BlankScreenSec*sampling_freq:end),'k|',5) %,BlankScreenSec*sampling_freq:end))
                xtickc=int32(linspace(1*10^-3*sampling_freq,length(Data.SortedRasters{AC(c)}{d}(:,BlankScreenSec*sampling_freq:end)),3)); % -BlankScreenSec*sampling_freq
                names= round(((double(xtickc)/sampling_freq)-10^-3),1)*10^3;
                set(gca, 'XTick',  xtickc, 'XTickLabel', names)
                axis square
                %xlim([4*sampling_freq 5*sampling_freq])
                title([num2str(Orientations(i)),char(176)]);
                
                %ylabel(['CPD: ',num2str(CPDs(p))],'FontSize',2);
                PSTHBineSize = 0.15; % binsize for PSTH in Seconds.
                %[PSTH{c}{d}{p}(i,:),binsize_sec]=Build_psth5(Data.SortedRasters{AC(c)}{d}(OrientationStimNums{i}(:,1),:),sampling_freq,PSTHBineSize);
                CountWindow = [2 22];
                [PSTH{c}{d}{p}(i,:),binsize_sec,smoothed_y{c}{d}{p}(i,:),SpikeCount{c}{d}{p}(i,:),Label]=Build_psth3(Data.SortedRasters{AC(c)}{d}(OrientationStimNums{i}(:,1),:),sampling_freq,CountWindow);
                %PSTHBinsForResponse = floor(StimTimeForPSTHCalc/binsize_sec);
                %ResponsePerOrientation{c}{d}(p,i) = max(PSTH{c}{d}{p}(i,(BlankScreenSec/binsize_sec):end)); % Calculate the max response per orientation over entire stimulus presentation in spikes/sec.
                %ResponsePerOrientation{c}{d}(p,i) = sum(PSTH{c}{d}{p}(i,(BlankScreenSec/binsize_sec):end))/(round(mean(diff(stimulus_times{RC})),1)-BlankScreenSec); % Calc Total Rate over Stim period. 
                CountPerOrientation{c}{d}(p,i) = sum(PSTH{c}{d}{p}(i,(BlankScreenSec/binsize_sec):end));
                count = count+1;
            end
   % Sponteneous Activity 
    for i=1:Data.dim{AC(c)}
    Data.SortedSpikesIdx{AC(c)}{i} = Aligned_idx(Data.ClusterIdx{AC(c)}==i);
    [Data.Spon{AC(c)}{i}] = OrientationSponCalc(Data.SortedSpikesIdx{AC(c)}{i},stimulus_times{AC(c)},sampling_freq,BlankScreenSec);
    end

            prompt = {'Count / Rate? (1=Count | 0=Rate)'};
            dlgtitle = 'Input';
            dims = [1 35];
            definput1 = {'1'};
            answer1 = inputdlg(prompt,dlgtitle,dims,definput1);
            if answer1{1} == '0'
            figure(OrResponsePerCPDFig{c}{d})
            subplot(1,length(CPDs),p)
            plot([0:30:330],ResponsePerOrientation{c}{d}(p,:),'b*');
            hold on
            ResponsePerCPD(c,d) = mean(ResponsePerOrientation{c}{d}(p,:));
            plot([0:30:330],ones(length([0:30:330]),1)*Data.Spon{AC(c)}{d},'r--')
            ylim([0 max(ResponsePerOrientation{c}{d}(p,:))+10]);
            xlabel('Orientation [Deg]');
            xticks([0:30:330]);
            xlim([-10 340]);
            ylabel('Firing Rate [Spikes/s]');
            %title(['CPD: ',num2str(CPDs(p))]);
            else
            figure(OrResponsePerCPDFig{c}{d})
            subplot(1,length(CPDs),p)
            plot([0:30:330],CountPerOrientation{c}{d}(p,:),'b*');
            hold on
            ResponsePerCPD(c,d) = mean(CountPerOrientation{c}{d}(p,:));
            plot([0:30:330],ones(length([0:30:330]),1)*Data.Spon{AC(c)}{d},'r--')
            ylim([0 max(CountPerOrientation{c}{d}(p,:))+1]);
            xlabel('Orientation [Deg]');
            xticks([0:30:330]);
            xlim([-10 340]);
            ylabel('Spike Count');
            end
            end
    end
Count = Count+1;
end
 %% OSI calculation
% for c = 1:length(AC)
%     for i=1:Data.dim{AC(c)}
%         for k=1:size(ResponsePerOrientation{c}{i},1)
%             MaxPerCPD(k) = max(ResponsePerOrientation{c}{i}(k,:));
%         end
%         [m,j] = max(MaxPerCPD);
%         %j = 2;
%         [Rpref,RprefIdx] = max(ResponsePerOrientation{c}{i}(j,:));
%         RprefDeg = (RprefIdx-1)*30;
%         if RprefDeg < 270
%             RorthDeg = RprefDeg+rad2deg(pi/2);
%             Rorthidx = RorthDeg/30+1;
%             Rorth = ResponsePerOrientation{c}{i}(j,Rorthidx);
%         else
%             RorthDeg = RprefDeg-rad2deg(pi/2);
%             Rorthidx = RorthDeg/30+1;
%             Rorth = ResponsePerOrientation{c}{i}(j,Rorthidx);
%         end
%         OSI{c}{i} = (Rpref-Rorth)/(Rpref+Rorth);
%     end
% end

