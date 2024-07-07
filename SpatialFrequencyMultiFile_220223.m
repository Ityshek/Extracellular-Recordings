%% Load the data files and Concat
clear all; close all;


AC = [str2num(cell2mat(inputdlg('Insert the Numbers of the Recorded Channels')))];
[raw_data, fs,stim_Data,stim_sampling_rate,Begin_record,stimulus_times,Stim_indx,CPDs,pathname] =load_data_ConcateMultiUnit(AC);

%% Plot Raw Data
std_factor = -4;        
FigRawData = figure();
for c = 1:length(AC)
    % Find Spikes
    [Data.Spike_ind,waveforms,thresh]=find_spikesANDwaves(AC,raw_data,fs,std_factor,Stim_indx);
    [aligned_Spikes,Data.Average_Spike,Aligned_idx]=Align_spikesRF(AC,waveforms,fs,std_factor,Data.Spike_ind);
    % Calc SNR
    [Data.SNR{AC(c)}] = SNRCalc(waveforms{AC(c)},raw_data{AC(c)}(1:Stim_indx{AC(c)}(1)),thresh(AC(c)),Data.Spike_ind{AC(c)});
    
    
    % Plot Raw Data + Spikes
    t=[0:length(raw_data{AC(c)})-1]/fs;
    threshold = thresh(AC(c))*ones(1,length(t));
    x = nan(1,length(t));
    x(Data.Spike_ind{AC(c)})=raw_data{AC(c)}(Data.Spike_ind{AC(c)});
    stim=nan(length(raw_data{AC(c)}),1); stim(Stim_indx{AC(c)})=90;
    figure(FigRawData);
    subplot(length(AC),1,c)
    plot(t,raw_data{AC(c)},'b',t,threshold,'r',t,x,'*k',t,stim,'vg')
    % figure settings
    xlabel('Time [Sec]','FontSize',20)
    ylabel('Amplitude [\muV]','FontSize',20)
    title(['Channel ',num2str(AC(c))])
    %legend('Raw Signal',['Threshold (SNR = ',num2str(Data.SNR{c}),')'],'Detected Spikes','Trigger');
    ylim([-200 200])
    xlim([0 max(t)])
    hold on
end
%% PCA & Clustering
%close all;
ActiveChannels = cell2mat(inputdlg('Select Channels for Further Analysis'));
ActiveChannels = str2num(ActiveChannels);
FigPlotNum = length(ActiveChannels);
AvgSpkFig = figure;  ClusterResultsFig = figure;
col=['r','g','m','c','y','k'];
for c=1:length(ActiveChannels)

    % Average Spike
    [Data.AlignedSpikes{ActiveChannels(c)},Data.AverageSpike{ActiveChannels(c)},Aligned_idx]=Align_spikes4(waveforms{ActiveChannels(c)},fs,std_factor, Data.Spike_ind{ActiveChannels(c)});
    t_spike=((0:length(Data.AverageSpike{ActiveChannels(c)})-1)/fs)*10^3;
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
    ClustEvalCH = evalclusters(Data.AlignedSpikes{ActiveChannels(c)},'kmeans','CalinskiHarabasz','KList',[1:2]);
    ClustEvalG = evalclusters(Data.AlignedSpikes{ActiveChannels(c)},'kmeans','gap','KList',[1:2]);
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
        t_sort=((0:size(Data.SortedSpikes{1,ActiveChannels(c)}{i},1)-1)/fs)*10^3;
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
        ylim([-250 250]); xlim([0.5 3]);
    end



    %  raw data with sorted spikes
    %         stim=nan(length(raw_data),1);
    %         stim(stimulus_indexes)=90;
    figure();
    subplot(FigPlotNum,FigPlotNum,c)
    %       threshold = Data.thresh(c)*ones(1,length(t));
    plot(t,raw_data{AC(c)},'b',t,stim,'vg',t,threshold,'r')

    % figure settings
    xlabel('Time[Sec]','FontSize',20)
    ylabel('Amplitude[\muV]','FontSize',20)
    title(['Channel ',num2str(ActiveChannels(c))])
    ylim([-400 400])
    xlim([0 max(t)])
    hold on

    for i=1:Data.dim{ActiveChannels(c)}
        Data.ClusteredspikeIdx{c}{i} = nan(1,length(Aligned_idx));
        Data.ClusteredVal{i} = nan(1,length(Aligned_idx));
        for k = 1:length(Aligned_idx)
            if  Data.ClusterIdx{ActiveChannels(c)}(k) == i
                Data.ClusteredspikeIdx{c}{i}(k) = Aligned_idx(k);
                Data.ClusteredVal{i}(k) = raw_data{AC(c)}(Aligned_idx(k));
            end
        end
        y{i} = nan(1,length(t));
        y{i}(rmmissing(Data.ClusteredspikeIdx{c}{i})) = rmmissing(Data.ClusteredVal{i});
        plot(t,y{i},['*';col(i)]);
        hold on
    end
    legend('Raw Signal','Trigger','Threshold','Detected Spikes Cluster 1','Detected Spikes Cluster 2')
    %   TIH
%     TIHFig = figure;
%     for i=1:Data.dim{ActiveChannels(c)}
%         Data.ISI{ActiveChannels(c)}{i} = (diff(rmmissing(Data.ClusteredspikeIdx{i}))/fs)*10^3; % Save the ISIs in ms for every cluster
%         figure(TIHFig)
%         subplot(2,2,i)
%         edges = [0:0.5:100];
%         histogram(Data.ISI{ActiveChannels(c)}{i}...
%         ,edges,'Normalization','pdf');
%         xlabel('ISI [ms]')
%         ylabel('Probability')
%         ylim([0 0.3])
%     end
   %CorrPlot = CorrFunc(t,Data.ClusteredspikeIdx,Data.dim{ActiveChannels(c)},sampling_freq);
end
%% Select Relevant Clusters
Data.Clusters = [str2num(cell2mat(inputdlg('Insert the Number of the Relevant Clusters')))]; % For multi file analysis, insert cluster numbers according of the unit's order from previous files.
%% Raster & PSTH & Curve
NumCPD = length(CPDs);
Nreps = 101;
BlankScreenSec = 1;
CountWindow = [102 152]; %win size in ms divided by 10
OffCountWindow = CountWindow-100;
for c=1:length(ActiveChannels)
FigCPDCurve = figure();
    for d=1:Data.Clusters
    CPDVec = [];
        % Raster Plots for all CPDs
FigRaster = figure(); FigPSTH = figure(); 
count = 1;
for IdxCPD = 1:Nreps:NumCPD*(Nreps)
    [Data.SortedRasters{AC(c),d}{count}] = build_rast_sort4(Data.ClusterIdx{AC(c)},Stim_indx{AC(c)}(IdxCPD:IdxCPD+Nreps-1)...
        ,fs,Data.dim{AC(c)},Aligned_idx,stimulus_times{AC(c)}(IdxCPD:IdxCPD+Nreps-1));
    figure(FigRaster);
    subplot(1,floor(NumCPD),count)
    spy(Data.SortedRasters{AC(c),d}{count}{1},'k|',5) %,BlankScreenSec*sampling_freq:end))
    xtickc=int32(linspace(1*10^-3*fs,length(Data.SortedRasters{AC(c),d}{count}{1}(:,:)),5)); % -BlankScreenSec*sampling_freq
    names= round(((double(xtickc)/fs)-10^-3),1)*10^3;
    set(gca, 'XTick',  xtickc, 'XTickLabel', names)
    axis square
    title(['CPD: ', num2str(CPDs(count))])
    % PSTH all CPDs
    
    [Data.PSTH{c}{d}{count},binsize_sec,Data.smoothedPSTH{c}{d}{count},Data.SpikeCount{c}{d}(count),Label]=Build_psth3(Data.SortedRasters{AC(c),d}{count}{1},fs,CountWindow);
    t_pst=1000*linspace(-10*10^-3,size(Data.PSTH{c}{d}{count},2)*binsize_sec-10*10^-3,...
        length(Data.PSTH{c}{d}{count}));
    
    % Off response Calc
    temp = full(Data.SortedRasters{AC(c),d}{count}{1}(:,OffCountWindow(1)*fs*0.01:OffCountWindow(2)*fs*0.01));
    Data.OffCount{c}{d}(count) = sum(sum(temp))/size(temp,1);
    
    % Plot
    figure(FigPSTH);
    subplot(1,floor(NumCPD),count)
    bar(t_pst,Data.PSTH{c}{d}{count})
    hold on
    plot(t_pst,Data.smoothedPSTH{c}{d}{count},'linewidth',2)
    xlabel('Time [ms]')%,'FontSize',20)
    ylabel(Label)%,'FontSize',20)
    ylim([0 0.8])
    title(['CPD: ', num2str(CPDs(count))])
    CPDVec = [CPDVec;CPDs(count)];
    
    % Calc Baseline Activity Every CPD
    if IdxCPD > 1
    SponSpikes = intersect(find(Data.ClusteredspikeIdx{c}{d} > Stim_indx{ActiveChannels(c)}(IdxCPD-1))...
    ,find(Data.ClusteredspikeIdx{c}{d} < Stim_indx{ActiveChannels(c)}(IdxCPD)));
    LastSponSpike = max(SponSpikes);

    SponSpikeTrain = zeros(1,int32(Data.ClusteredspikeIdx{c}{d}(LastSponSpike)-Stim_indx{ActiveChannels(c)}(IdxCPD-1)));
    SponSpikeTrain(int32(Data.ClusteredspikeIdx{c}{d}(SponSpikes))) = 1; % Find all spikes of the current cluster before stimulus onset
   
    BinSizeSpon = (CountWindow(2) - CountWindow(1))*440; countSpon = 1;
    for k=BinSizeSpon:BinSizeSpon:size(SponSpikeTrain,2) - BinSizeSpon
        SponBinned(countSpon) = sum(SponSpikeTrain(k-BinSizeSpon+1:k+BinSizeSpon+1));  % Calc num of Spon spikes in 200ms window (should be same window as Spike Count Var).
        countSpon = countSpon+1; 
    end
    Data.Spon{c}{d}(count) = mean(SponBinned);
    Data.SponStd{c}{d}(count) = std(SponBinned);
    else
        SponSpikes = find(Data.ClusteredspikeIdx{c}{d} > Stim_indx{ActiveChannels(c)}(1));
            LastSponSpike = max(SponSpikes);
    SponSpikeTrain = zeros(1,int32(Stim_indx{ActiveChannels(c)}(1)));
    SponSpikeTrain(int32(Data.ClusteredspikeIdx{c}{d}(SponSpikes))) = 1; % Find all spikes of the current cluster before stimulus onset
   
    BinSizeSpon = (CountWindow(2) - CountWindow(1))*440; countSpon = 1;
    for k=BinSizeSpon:BinSizeSpon:size(SponSpikeTrain,2) - BinSizeSpon
        SponBinned(countSpon) = sum(SponSpikeTrain(k-BinSizeSpon+1:k+BinSizeSpon+1));  % Calc num of Spon spikes in Count window.
        countSpon = countSpon+1; 
    end
    Data.Spon{c}{d}(count) = mean(SponBinned);
    Data.SponStd{c}{d}(count) = std(SponBinned); 
    end
    clear SponSpikeTrain LastSponSpike SponSpikes SponBinned
    count = count+1;
end

% Plot Curve
figure(FigCPDCurve);
SponFlag = inputdlg('Subtracte Baseline Activity? (1 = yes | 0 = No)');
Data.response_magnitude = mean(Data.SpikeCount{c}{1}/Data.Spon{c}{d});
if SponFlag{1} == '1'
    Data.CountNoSpon{c}{d} = Data.SpikeCount{c}{1}-Data.Spon{c}{d};
    for k=1:length(Data.CountNoSpon{c}{d})
        if Data.CountNoSpon{c}{d}(k) < 0
            Data.CountNoSpon{c}{d}(k) = 0;
        end
    end
plot(CPDVec,Data.CountNoSpon{c}{d})
else
plot(CPDVec,Data.SpikeCount{c}{d})
end
xlabel('CPD'); ylabel(['Spike Count [',num2str((CountWindow(2)-CountWindow(1))*10),'ms]'])
end
end
%% Calc % Over Spon
% Second Variable in {} should be the relevant cluster
Data.spon_per = max(Data.SpikeCount{1}{Data.Clusters})/mean(Data.Spon{1}{Data.Clusters})*100;
%% Save the Plots
% Change "SavePath" to the desiered folder
[a,b] = regexp(pathname,'Data\\\S{10,10}');
date = pathname(a+5:b);
[a,b] = regexp(fname,'Loc\d*'); Loc = fathname(a:b);
SavePath = ('C:\Users\....\Figures\');
save([SavePath,'SpatialFrequencyPlot_',date,Loc,'.fig'],"FigCPDCurve");
%% Save the Data Variable
% Change "SavePath" to the desiered folder
SavePath ='/Users/shirahasky/Desktop/practicum/Data/';
% stimfreq = num2str(round(1/mean(diff(stimulus_times)))); 
% [a,b] = regexp(fname,'_\d*ms'); Data.StimDur = fname(a+1:b-2); 

[a,b] = regexp(pathname,'/[1234567890.]*/'); date = pathname(a+1:b-1);
[a,b] = regexp(pathname,'/[1234567890Loc]*/'); Location = pathname(a+1:b-1);
Ans = mat2str(CPDs);
Ans = strrep(Ans,' ','_');
Data.CPDs = Ans(2:end-1);
save([SavePath,'SpatialFrequency_',Ans(2:end-1),'CPD',date,Location,'Channel',num2str(AC),'.mat'],"Data");