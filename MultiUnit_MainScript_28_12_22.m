%% Load & Plot Basic Data
close all;
clearvars -except Data;
RawFig = figure;
%ChannelPosition16 = [1,9,5,6,14,13,2,10,15,7,12,11,3,4,16,8];
%OpenChannels = ['3,7,11,15,19,23,27,31,2,6,10,14,18,22,26,30,4,8,12,16,20,24,28,32,1,5,9,13,17,21,25,29'];
%OpenChannels = ['5,9,13,17,21,25,29'];
%NumChannels = length(split(OpenChannels,','));
std_Factor=-4.5;
prompt = {'Select Projection System (Screen = 1, Micron =2)', 'Number of Recorded Channels:','All Desired Channels:', 'N rows:','N Column:', 'N Reps:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'2','4','1,4,9,13','2','2','50'};
%definput = {'2','31',['1,2,3,4,5,6,7,8,9,10,11,12,13,15,16...' ...
%   '17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32'],'4','8','50'};
%definput = {'2',num2str(NumChannels),num2str(OpenChannels),'3','3','80'}
answer = inputdlg(prompt,dlgtitle,dims,definput);
FigPlotNum = str2double(answer{2}); count = str2double(answer{3}); nRows = str2num(answer{4}); nColumn = str2num(answer{5});
definptOut = {'150'};
if answer{1} == '2'
    prompt = {'Insert Projection Type (1=NIR | 2=Vis | 3=HDMI)','Is the trigger signal distorted? (1=YES, 0=NO)'...
        ,'Insert Stim Type: [Alter = 1,CFF = 2, Intensity = -]'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput1 = {'2','0','2'};
    answer1 = inputdlg(prompt,dlgtitle,dims,definput1);
end
pathname = []; SignalFiles = []; Dir = [];

if answer1{1} ~= '3'
    [fname,pathname]=uigetfile('*.mat','Choose data file');
end

chan = split(answer{3},',');
if answer1{1} ~= '3'
    RastFig = figure(); PsthFig = figure();
end
flag = []; % indicates if a channel was loaded already
for c = 1:FigPlotNum
    %RC = [str2double(cell2mat(inputdlg('Insert the Number of the Recorded Channel',dlgtitle,dims,chan(c))))];
    RC =str2double(chan(c));
    % SELECT ONE OF THE TWO ROWS BELOW ACCORDING TO THE STIMULATION SYSTEM
    if answer{1} == '1'
        [raw_data{c}, sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag] =load_data_Screen(RC,fname,pathname); % data loader for Screen stimulation
        for i = 1:length(raw_data{1}); if abs(raw_data{1}(i)) > 500; raw_data{1}(i) = mean(raw_data{1}); end; end
        [stimulus_times,stimulus_indexes]=find_stim(stim_Data,stim_sampling_rate,sampling_freq,Begin_record); % run only for screen stimulus
        VisFlag = cell(1);
    else
        if answer1{1} == '3'
            [raw_data{c},sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag,stimulus_times, stimulus_indexes,VisFlag,fname,pathname,SignalFiles,Dir] = load_data_MicronHDMI(RC,answer1,flag,pathname,SignalFiles,Dir); % data loader for Micron HDMI stimulation
            fname = fname{1};
            flag = 1;
        else
            [raw_data{c},sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag,stimulus_times, stimulus_indexes,VisFlag] = load_data_Micron(RC,fname,pathname,answer1); % data loader for Micron stimulation
        end
    end

    if sum(ismember('Reps',fname)) == 4
        [startIndex,endIndex] = regexp(fname,'_[0123456789._]*Reps');
        NumReps = str2num(fname(startIndex+1:endIndex-4));
    else
        NumReps = str2num(answer{6});
    end
    if ~channelflag


        %         [startIndex,endIndex] = regexp(fname,'_\d*ms');
        %         if startIndex ~= 0
        %             StimDuration = str2double(fname(startIndex+1:endIndex-2));
        %         else
        %             StimDuration = 0;
        %         end
        %stimulus_indexes = stimulus_indexes+round(StimDuration/1000*stim_sampling_rate);
        %stimulus_times = stimulus_times - (1-StimDuration*10^-3);
        t=[0:length(raw_data{c})-1]/sampling_freq;
        Data.thresh(c)=std_Factor*nanstd(double(raw_data{c}));
        %ZoomRaw = [60 62];
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
        ylim([-150 150])
        box off
        set(gca,'color','none')
        % Build raster%
        prompt = {'Insert Upper Thershold Value for spike outliers:'};
        dlgtitle = 'Input'; dims = [1 35]; Check=[]; PreStimTime = 200;
        %outlier = str2num(str2mat(inputdlg(prompt,dlgtitle,dims,definptOut)));
        outlier = str2num(str2mat(definptOut));
        % Detect Spikes, remove Outlier and build Raster 
        [Rast,Spike,Av_spike,indx_spike,ind_rast,spike_stim,spike_times]=build_rastRef(stimulus_indexes,stimulus_times,raw_data{c},sampling_freq,Data.thresh(c),outlier,VisFlag,Check,PreStimTime);
        Data.SpikeTimes{c} = spike_times;
        Data.IndxSpike{c} =  indx_spike;
        Data.Spike{c} = Spike;
        if answer1{1} == '3' % Check for Alternating Stim Condition
            [a,b] = regexp(fname,'_\d*OR'); NumOrs = str2num(fname(a+1:b-2)); % Extract number of orientations from file name
            RastFig = figure;
            tiledlayout('flow');
            count = 1;
            for r = 1:NumReps:NumOrs*NumReps
                AltRepsSplit(count,:) = round(linspace(r,r+NumReps,NumReps));
                count = count+1;
            end
            for r = 1:NumOrs % Split Rast to Orientations
                RastAlt{r} = Rast(AltRepsSplit(r,:),:);
                nexttile
                spy(RastAlt{r});
                % Raster figure settings
                xtickc = linspace(1,round(median(diff(stimulus_times(2:end))),2)*sampling_freq,11);
                names= round(((xtickc/sampling_freq)),2)*10^3-PreStimTime;
                set(gca, 'XTick',  xtickc, 'XTickLabel', names)
                xlim([-440 max(xtickc)]);
                axis square
                %xlabel('Time [ms]','FontSize',12)
                %ylabel('Stimulus Repetition','FontSize',12)
                title(['Channel ',num2str(RC),' ',num2str(45*(r-1)),'Deg'])
                count = count+1;
            end
                
        else
            figure(RastFig);
            subplot(nRows,nColumn,c)
            spy(Rast);
            %x=(set(gca,'XTickLabel',(linspace(-PreStimTime*10^-3,mean(diff(stimulus_times))-PreStimTime*10^-3,5))));
            %             xtickc=linspace(-PreStimTime*10^-3*sampling_freq,round(median(diff(stimulus_times(2:end))),2)*sampling_freq...
            %                 -PreStimTime*10^-3*sampling_freq,5);
            xtickc = linspace(1,round(median(diff(stimulus_times(2:end))),2)*sampling_freq,11);
            names= round(((xtickc/sampling_freq)),2)*10^3-PreStimTime;

            % figure settings
            set(gca, 'XTick',  xtickc, 'XTickLabel', names)
            xlim([0 max(xtickc)]);
            axis square
            %xlabel('Time [ms]','FontSize',12)
            %ylabel('Stimulus Repetition','FontSize',12)
            title(['Channel ',num2str(RC)])
        end



        % Build PSTH%
        prompt = {'Choose Spike Calc: (1 = Hz | 2 = Count)','Choose Window Size [ms]:'};
        definpt = {'2','10'}; dlgtitle = 'Input'; dims = [1 35];
        %PSTHSettings = inputdlg(prompt,dlgtitle,dims,definpt);
        PSTHSettings = definpt;
        if answer1{1} == '3' % Check for Alternating Stim Condition
            CountWindow = [2 22];
            PsthFig = figure;
            tiledlayout('flow');
            count = 1;
            for r = 1:length(RastAlt) % Loop orientations
                %ax = axes();
                [PsthBinned,binsize_sec,smoothed_Psth,Label]=Build_psth3(RastAlt{r},sampling_freq,CountWindow,PSTHSettings);
                %StimTime = round(mean(diff(stimulus_times)),2);
                t_pst=1000*linspace(-PreStimTime*10^-3,size(PsthBinned,2)*binsize_sec-PreStimTime*10^-3,length(PsthBinned));
                nexttile
                b = bar(t_pst,PsthBinned,DisplayName=[num2str(binsize_sec*1000) 'ms Bins']);
                hold on
                plot(t_pst,smoothed_Psth,'linewidth',2,DisplayName='5 Bin smoothning')
                % figure settings
                %xticklabels(linspace(0,StimTime*1000,5))
                %xticklabels([-10])
                %xlabel('Time [ms]','FontSize',12)
                %ylabel(Label,'FontSize',12)
                title(['Channel ',num2str(RC),' ',num2str(180/NumOrs*(r-1)),'Deg'])
                ylim([0 1]);
                xlim([-PreStimTime max(t_pst)]);
                axis square; box off;
                set(gca,'color','none','FontSize',15)
                %       ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
                %       ax.Box = 'off'; ax.Color = "none";
                
            end
        else
            CountWindow = [2 22];
            figure(PsthFig);
            ax = axes();
            subplot(nRows,nColumn,c)
            [PsthBinned,binsize_sec,smoothed_Psth,Label]=Build_psth3(Rast,sampling_freq,CountWindow,PSTHSettings);
            Data.PSTH{c} = PsthBinned;
            StimTime = round(mean(diff(stimulus_times)),2);
            t_pst=1000*linspace(-PreStimTime*10^-3,size(PsthBinned,2)*binsize_sec-PreStimTime*10^-3,length(PsthBinned));
            b = bar(t_pst,PsthBinned,DisplayName=[num2str(binsize_sec*1000) 'ms Bins']);
            hold on
            plot(t_pst,smoothed_Psth,'linewidth',2,DisplayName='5 Bin smoothning')
            % figure settings
            %xticklabels(linspace(0,StimTime*1000,5))
            %xticklabels([-10])
            %xlabel('Time [ms]','FontSize',12)
            %ylabel(Label,'FontSize',12)
            title(['Channel ',num2str(RC)])
            ylim([0 1.5]);
            xlim([-PreStimTime max(t_pst)]);
            axis square; box off;
            set(gca,'color','none','FontSize',15)
            %       ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
        end
    end
    count = count+1;
    figure(RawFig);
    x = nan(1,length(t));
    x(indx_spike)=raw_data{c}(indx_spike);
    hold on
    plot(t,x,'*k')
    %xlabel('Time[Sec]','FontSize',12)
    %ylabel('Amplitude[\muV]','FontSize',12)
    %% Calc Noise
    Sections = 25;
    for n=1:Sections
        NoiseMean(n) = length(find(spike_times<((stimulus_times(1)/Sections)*n)));
        if n>1
            NoiseMean(n) = NoiseMean(n)-sum(NoiseMean(1:n-1));
        end
    end
    NoiseMean = (NoiseMean/(stimulus_times(1)/Sections))*binsize_sec;
    NoiseLevel = mean(NoiseMean)+2*std(NoiseMean);
    % Plot Noise Level
    % figure(PsthFig);
    % plot(t_pst,NoiseLevel*ones(1,length(t_pst)),'k-',LineWidth=2,DisplayName='Noise Level + 2std')
    % legend();
    % Calculate SNR
    %     SponIdxSpikes = indx_spike(find(indx_spike<stimulus_indexes(2)));
    %     SponSpike = Spike(1:length(SponIdxSpikes));
    %     [Data.SNR{c}] = SNRCalc(SponSpike,raw_data{c}(1:stimulus_indexes(2)),Data.thresh(c),SponIdxSpikes);
end
    figure(RawFig);suplabel('Time [Sec]','x'); suplabel('Amplitude[\muV]','y');
    figure(PsthFig); suplabel('Time [ms]','x'); suplabel(Label,'y');
    figure(RastFig); suplabel('Time [ms]','x'); suplabel('Stimulus Repetition','y');
%figure(RawFig);
%xlabel('Time[Sec]','FontSize',12)
%ylabel('Amplitude[\muV]','FontSize',12)
% legend('Raw Signal','Trigger',['Threshold (SNR = ',num2str(Data.SNR{c}),')'],'Detected Spikes');

%% Plot Stim Channel + Trigger Points
figure();plot(stim_Data)
hold on; plot(stimulus_times*2750,ones(1,size(stimulus_times,2))*(-1950),'*')
title('Stimulus Channel')
figure();plot(round(diff(stimulus_times),4));
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
    %Data.dim{ActiveChannels(c)}= 3;

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

end
%% Refractory Period Check
% TIH
TIHFig = figure;
figure(TIHFig)
tiledlayout('flow')
for i=1:Data.dim{ActiveChannels(c)}
    Data.ISI{ActiveChannels(c)}{i} = (diff(rmmissing(Data.ClusteredspikeIdx{i}))/sampling_freq)*10^3; % Save the ISIs in ms for every cluster
    nexttile
    edges = [0:1:50];
    histogram(Data.ISI{ActiveChannels(c)}{i}...
        ,edges,'Normalization','probability');
    xlabel('ISI [ms]')
    ylabel('Probability')
    ylim([0 0.3])
end
%CorrPlot = CorrFunc(t,Data.ClusteredspikeIdx,Data.dim{ActiveChannels(c)},sampling_freq);
%% Sorted Plots
clear  PsthBinned
for c = 1:length(ActiveChannels)
    [Rast_sort,ind_rast_sort]=build_rast_sort2(Data.ClusterIdx{ActiveChannels(c)},stimulus_indexes,sampling_freq,ind_rast,Data.dim{ActiveChannels(c)},stimulus_times);
    for i=1:Data.dim{ActiveChannels(c)}
        if answer1{1} == '3' % Check for Alternating Stim Condition
            [a,b] = regexp(fname,'_\d*OR'); NumOrs = str2num(fname(a+1:b-2)); % Extract number of orientations from file name
            RastSortFig = figure;
            tiledlayout('flow');
            count = 1;
            for r = 1:NumReps:NumOrs*NumReps % Loop orientations
                RastSortAlt{i}{count} =Rast_sort{i}(r:r+NumReps-1,:); % Extract relevant repetitions from Raster
                nexttile
                spy(RastSortAlt{i}{count})

                x=(set(gca,'XTickLabel',(linspace(-PreStimTime*10^-3,mean(diff(stimulus_times))-PreStimTime*10^-3,5))));
                xtickc=linspace(-PreStimTime*10^-3*sampling_freq,round(median(diff(stimulus_times(2:end))),2)*sampling_freq...
                    -PreStimTime*10^-3*sampling_freq,5);
                names= round(((xtickc/sampling_freq)-PreStimTime^-3),2)*10^3;

                % figure settings
                set(gca, 'XTick',  xtickc, 'XTickLabel', names)
                xlim([-PreStimTime max(xtickc)]);
                axis square
                xlabel('Time [ms]','FontSize',15)
                ylabel('Stimulus Repetition','FontSize',15)
                title(['Unit ',num2str(i),' ',num2str(45*(count-1)),'Deg'])
                count = count+1;
            end
        else
            sorted_rast(i)=figure;
            spy(Rast_sort{i})
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
    end




    % PSTH Sorted
    for i=1:Data.dim{ActiveChannels(c)}
        prompt = {'Choose Spike Calc: (1 = Hz | 2 = Count)','Choose Bin Size [ms]:'};
        definpt = {'2','20'}; dlgtitle = 'Input'; dims = [1 35];
        PSTHSettings = inputdlg(prompt,dlgtitle,dims,definpt);
        window = [20/str2num(PSTHSettings{2}) 320/str2num(PSTHSettings{2})];
        if PSTHSettings{1} == '1' % Choose Y axis label acording to calc method
            Label = 'Firing Rate [Spikes/sec]';
        else
            Label = ['Spike Count [',num2str(PSTHSettings{2}),'ms]'];
        end
        if answer1{1} == '3' % Check for Alternating Stim Condition
            PsthFig{i} = figure;
            tiledlayout('flow');
            count = 1;
            for r = 1:NumReps:NumOrs*NumReps % Loop orientations
                %ax = axes();
                [PsthBinned{i}{count},binsize_sec,Psth_sort{i}{count}]=Build_psth4(RastSortAlt{i}{count},sampling_freq,window);
                %StimTime = round(mean(diff(stimulus_times)),2);
                t_pst=1000*linspace(-PreStimTime*10^-3,size(PsthBinned{i}{count},2)*binsize_sec-PreStimTime*10^-3,length(PsthBinned{i}{count}));
                nexttile
                b = bar(t_pst,PsthBinned{i}{count},DisplayName=[num2str(binsize_sec*1000) 'ms Bins']);
                hold on
                plot(t_pst,Psth_sort{i}{count},'linewidth',2,DisplayName='5 Bin smoothning')
                % figure settings
                %xticklabels(linspace(0,StimTime*1000,5))
                %xticklabels([-10])
                xlabel('Time [ms]','FontSize',15)
                ylabel(Label,'FontSize',15)
                title(['Unit ',num2str(i),' ',num2str(45*(count-1)),'Deg'])
                ylim([0 50]);
                xlim([-PreStimTime max(t_pst)]);
                axis square; box off;
                set(gca,'color','none','FontSize',15)
                %       ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
                %       ax.Box = 'off'; ax.Color = "none";
                count = count+1;
            end
        else
            [PsthBinned(i,:),binsize_sec,Psth_sort(i,:)]=Build_psth4(Rast_sort{i},sampling_freq,window,PSTHSettings);
            t_pst=1000*linspace(-PreStimTime*10^-3,size(Psth_sort(i,:),2)*binsize_sec-PreStimTime*10^-3,length(Psth_sort(i,:)));
            sorted_PSTH(i)=figure;
            bar(t_pst,PsthBinned(i,:))
            hold on
            plot(t_pst,Psth_sort(i,:),'linewidth',2)
            xlabel('Time [ms]','FontSize',20)
            ylabel(Label,'FontSize',20)
            axis square; box off;
            set(gca,'color','none','FontSize',20)
            xlim([-PreStimTime max(t_pst)]); ylim([0 1.5]);
            title(['Sorted PSTH ', num2str(i)])
        end
    end

end


%% Select Relevant Clusters
Data.Clusters = [str2num(cell2mat(inputdlg('Insert the Number of the Relevant Clusters')))]; % For multi file analysis, insert cluster numbers according of the unit's order from previous files.
if ~isfield(Data,'Rast_sort')
    Data.Rast_sort = {}; Data.Psth_sort = {}; Data.ind_rast_sort ={}; Data.PsthBinned = {}; Data.WaveformMat = {};
    Data.SpikeIdx = {};
    if answer1{1} == '3' % Check for Alternating Stim Condition
        Data.WaveformMat{end+1} = Data.SortedSpikes{Data.Clusters};
        Data.SpikeIdx{end+1} = rmmissing(Data.ClusteredspikeIdx{Data.Clusters});
        Data.Rast_sort{end+1} = RastSortAlt{Data.Clusters};
        Data.Psth_sort{end+1} = Psth_sort{Data.Clusters};
        Data.ind_rast_sort{end+1} = ind_rast_sort{Data.Clusters};
        Data.PsthBinned{end+1} = PsthBinned{Data.Clusters};
    else
        Data.WaveformMat{end+1} = Data.SortedSpikes{Data.Clusters};
        Data.SpikeIdx{end+1} = rmmissing(Data.ClusteredspikeIdx{Data.Clusters});
        Data.Rast_sort{end+1} = Rast_sort{Data.Clusters};
        Data.Psth_sort{end+1} = Psth_sort(Data.Clusters,:);
        Data.ind_rast_sort{end+1} = ind_rast_sort{Data.Clusters};
        Data.PsthBinned{end+1} = PsthBinned(Data.Clusters,:);
    end
else
    if answer1{1} == '3' % Check for Alternating Stim Condition
        Data.WaveformMat{end+1} = Data.SortedSpikes{Data.Clusters};
        Data.SpikeIdx{end+1} = rmmissing(Data.ClusteredspikeIdx{Data.Clusters});
        Data.Rast_sort{end+1} = RastSortAlt{Data.Clusters};
        Data.Psth_sort{end+1} = Psth_sort{Data.Clusters};
        Data.ind_rast_sort{end+1} = ind_rast_sort{Data.Clusters};
        Data.PsthBinned{end+1} = PsthBinned{Data.Clusters};
    else
        Data.WaveformMat{end+1} = Data.SortedSpikes{Data.Clusters};
        Data.SpikeIdx{end+1} = rmmissing(Data.ClusteredspikeIdx{Data.Clusters});
        Data.Rast_sort{end+1} = Rast_sort{Data.Clusters};
        Data.Psth_sort{end+1} = Psth_sort(Data.Clusters,:);
        Data.ind_rast_sort{end+1} = ind_rast_sort{Data.Clusters};
        Data.PsthBinned{end+1} = PsthBinned(Data.Clusters,:);
    end
end
%% Sponteneous Activity calculated from recording prior to 1st trigger.
% Run next line Only for first trial of each experiment
prompt = {'Is First Trial? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0'};
answer2 = inputdlg(prompt,dlgtitle,dims,definput);
if answer2{1} == '1'
    Data.IsActive = {}; Data.ResponseWidth = {}; Data.Latency = {};
    Data.Spon = {}; Data.SponStd = {};
end
prompt = {'Choose Baseline Calc Method (1=Prior to 1st Trigger || 0 = Prior to each Trigger)',...
    'Choose Window Size for Response Calc',...
    'Choose Num of Bins for Activation Calc'};
dlgtitle = 'Input'; dims = [1 35];

% Window for Spike count Calcs:
%NFF = 200ms
%PFF = 100ms
%WindowSize = (window(2) - window(1))*binsize_sec;
definput = {'0','200','2'}; answer3 = inputdlg(prompt,dlgtitle,dims,definput);
SponFlag = answer3{1}; WindowSize = str2num(answer3{2});
Data.PSTHBinSize = binsize_sec; Data.WindowSize = WindowSize;
Data.ThreshBinSize = str2num(answer3{3});
% SponFlag = cell2mat(inputdlg('Choose Natural / Prosthetic (1 = N | 0 = P)'));
% if SponFlag == '0'

LastSponSpike = max(find(spike_times<stimulus_times(1)));
%if answer1{1} == '3' % Check for Alternating Stim Condition
% FirstOrs = linspace(1,NumReps*NumOrs-NumReps+1,NumOrs); % Find Idxs of first Reps for each Orientation
% FirstOrsIdx = stimulus_indexes(FirstOrs); % Find Stimulus Idxs for first reps.
% SponWindow = 20*sampling_freq; % Define time window to look for spon spikes (in seconds).
SponRelevantSpikes = rmmissing(Data.ClusteredspikeIdx{Data.Clusters}); %Extract the spike idxs of the relevant cluster
SponWindow = (window(2)-window(1))*0.01;
[Spon,SponStd] = SponCalc(Data.PsthBinned{end},Rast_sort{Data.Clusters},sampling_freq,stimulus_indexes,window,raw_data,SponRelevantSpikes,SponFlag,binsize_sec,WindowSize);
Data.Spon{end+1} = Spon; Data.SponStd{end+1} = SponStd;

%% Calc Activation Threshold
RelevantWindow = [(PreStimTime/(binsize_sec*1000))+1 ...
    ((PreStimTime/(binsize_sec*1000))+1)+(WindowSize/(binsize_sec*1000))]; % Create a window to look for response in PSTH Bins
% Data.IsActive{end+1} = zeros(1,length(Data.Rast_sort{end}));
% Data.ResponseWidth{end+1} = zeros(1,length(Data.Rast_sort{end}));
Data.IsActive{end+1} = []; Data.ResponseWidth{end+1} = [];
Data.Latency{end+1} = [];

if answer1{1} == '3' % Check for Alt Condition
    RelevantWindow2 = RelevantWindow+(1/binsize_sec);
    Data.IsActive{end+1} = zeros(1,length(Data.Rast_sort{end}));
    Data.ResponseWidth{end+1} = zeros(1,length(Data.Rast_sort{end}));
    Data.Latency{end+1} = zeros(1,length(Data.Rast_sort{end}));
    for r=1:length(Data.Rast_sort{end})
        RelevantBins2 = PsthBinned{Data.Clusters}{r}(RelevantWindow2(1):RelevantWindow2(2));
        A = find(RelevantBins2>=Spon+2*SponStd);
        if length(diff(find(A))) > Data.ThreshBinSize
            if find(diff(A) == 1)
                Data.IsActive{end}(r) = 1;
                Data.ResponseWidth{end}(r) = length(find(diff(A) == 1))*Data.PSTHBinSize*1000; % Calc Response Width in ms.
            end
        end
    end
else % Calc for Every Other Trial
    RelevantBins = Data.PsthBinned{end}(RelevantWindow(1):RelevantWindow(2)); % PSTH bins that are relevant for response
    A = find(RelevantBins>=Spon+2*SponStd);
    if length(diff(find(A)))+1 >= Data.ThreshBinSize
        if find(diff(A) == 1)
            Data.IsActive{end} = 1;
            Data.ResponseWidth{end} = (length(find(diff(A) == 1))+1)*Data.PSTHBinSize*1000; % Calc Response Width in ms.
            Data.Latency{end} = t_pst(A(1)+RelevantWindow(1)-1);
        end
    end
end
%% Intensity Response Calc
prompt = {'First or Last Trial? (1=1st;2=Last;0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0'};
answer4 = inputdlg(prompt,dlgtitle,dims,definput);

if answer4{1} == '1'
    Data.Intensity = [];
    Data.Response = {};
    Data.MaxResponse = {};
end
% Extract Intensity Values
if ~isnan(regexp(fname,'_[0123456789.]*CD')) % Screen Sim Intensity
    [Int1,Int2] = regexp(fname,'_[0123456789.]*CD');
    Data.Intensity = [Data.Intensity;str2num(fname(Int1+1:Int2-2))]; % Stimulus Intensity in CD/m^2.
elseif ~isnan(regexp(fname,'_[0123456789.]*nW')) % Micron Natural Intensity
    [Int1,Int2] = regexp(fname,'_[0123456789.]*nW'); StimType = 1;
    Data.Intensity = [Data.Intensity;str2num(fname(Int1+1:Int2-2))]; % Stimulus Intensity in nW/mm^2.
elseif ~isnan(regexp(fname,'FFNIR')) % Micron NIR Intensity
    if ~isnan(regexp(fname,'uW'))
        [a,b] = regexp(fname,'_[.0123456789]*uW'); StimType = 2;
    elseif ~isnan(regexp(fname,'Int'))
        [a,b] = regexp(fname,'_[.0123456789]*Int'); StimType = 2;
    end
    Data.Intensity = [Data.Intensity;str2num(fname(a+1:b-3))]; % Stimulus Intensity in uW/mm^2.
end
% Calc Responses
RastRelevantWindow = (RelevantWindow*binsize_sec - binsize_sec)*sampling_freq;
B=[];
for r = 1:size(Data.Rast_sort{end},1)
    B = [B,full(sum(Data.Rast_sort{end}(r,RastRelevantWindow(1):RastRelevantWindow(2))))];
end
Data.Response{end+1} = mean(B)*(1000/WindowSize); % Calc Mean Response in Spikes/Sec.
Data.MaxResponse{end+1} = max(Data.Psth_sort{end}(RelevantWindow(1):RelevantWindow(2))); % Calc Max Response in Spikes/Sec.
% Plot Relevant Graphs after all trials calculated
if answer4{1} == '2'
    % Plot the PSTH across all Stim Intensities
    Data.PsthsFig = figure();
    tiledlayout('flow')
    for p = 1:length(Data.Psth_sort)
        nexttile
        plot(t_pst,Data.Psth_sort{p},'b',LineWidth=2)
        hold on
        thresh = ones(1,length(t_pst))*(Data.Spon{p}+2*Data.SponStd{p});
        plot(t_pst,thresh,'r-',LineWidth=2)
        xlabel('Time[ms]');ylabel('Spike Rate[Hz]');
        xlim([t_pst(1) t_pst(end)+10]); ylim([0 100]);
        if StimType == 1
            title([num2str(Data.Intensity(p)),'nW/mm^2'])
        else
            title([num2str(Data.Intensity(p)),'mW/mm^2'])
        end
        axis square; box off;
        set(gca,'color','none','FontSize',15)
    end
    % Plot the Latency across all Stim Intensities
    Activevector = cellfun(@(x) ~isempty(x), Data.IsActive);
    Data.LatencyFig = figure();
    plot(Data.Intensity(find(Activevector == 1))',cell2mat(Data.Latency),'b',LineWidth=2)
    ylim([0 150]);
    ylabel('Response Latency [ms]');
    if StimType == 1
        Labelx =  ['Stimulus Intensity [nW/mm^2]']; xlim([10 650]);
    else
        Labelx =  ['Stimulus Intensity [mW/mm^2]']; xlim([0 2]);
    end
    xlabel(Labelx);
    axis square; box off;
    set(gca,'color','none','FontSize',15)
    % Plot Response across all stim Intensities
    Data.ResponseFig = figure();
    plot(Data.Intensity,cell2mat(Data.Response),LineWidth=2)
    hold on
    A = cell2mat(Data.Response); B = [find(Activevector == 0),max(find(Activevector == 0))+1];
    plot(Data.Intensity(B),A(B),'r-',LineWidth=2)
    xlabel(Labelx); ylabel('Spike Rate[Hz]');
    ylim([0 50]);
    axis square; box off;
    set(gca,'color','none','FontSize',15)
end
%% CFF Response
CFFCalcWindow = 1*sampling_freq; % Window to count the spikes in for response
SponToResponseFactor = CFFCalcWindow/(SponWindow*sampling_freq);
[a,b] = regexp(fname,'_[0123456789]*ms'); StimDuration = str2num(fname(a+1:b-2));
[Data] = CFFCalc(Data,CFFCalcWindow,Rast_sort{Data.Clusters},SponToResponseFactor,StimDuration,fname,sampling_freq);
%% Alternating Grid Different Orientations
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer2 = inputdlg(prompt,dlgtitle,dims,definput);

if answer2{1} == '1'
    Data.CPDs = [];  Data.Orientations = linspace(0,180-(180/NumOrs),NumOrs);
    Data.ResponseAlt = []; Data.Response = [];
end
[a,b] = regexp(fname,'_[.0123456789]*CPD');
Data.CPDs(end+1) = str2num(fname(a+1:b-3));
for r = 1:NumOrs
    [val,idx] = max(Data.PsthBinned{end}{r}(RelevantWindow(1):RelevantWindow(2))); % Find the idx of the peak in the relevant window
    Response(r) = mean(Data.PsthBinned{end}{r}(RelevantWindow(1)+idx-2:RelevantWindow(1)+idx)) - (Spon+2*SponStd); % Calc mean of 3 bins around the peak and reduce thresh value

    [val2,idx2] = max(Data.PsthBinned{end}{r}(RelevantWindow(1):RelevantWindow(2))); % Find the idx of the peak in the relevant window
    ResponseAlt(r) = mean(Data.PsthBinned{end}{r}(RelevantWindow2(1)+idx2-2:RelevantWindow2(1)+idx2)) - (Spon+2*SponStd); % Calc mean of 3 bins around the peak and reduce thresh value
end
Data.ResponseAlt = [Data.ResponseAlt; ResponseAlt]; Data.Response = [Data.Response; Response]; % Save Response to trials

if answer2{2} == '1'
    MaxResponse = max(max([Data.AltResponse,Data.Response]));
    if isnan(findstr('VIS',fname))
        Data.CPM = 107*Data.CPDs; % Calculated for NIR from visual inspection in Micron Camera.
    else
        Data.CPM = 48*Data.CPDs-0.3; % Calculated for VIS from visual inspection in Micron Camera.
    end
    Data.ResponseContourFig = figure();
    tiledlayout(2,1);
    ax1 = nexttile;
    contourf(Data.Orientations,Data.CPM,(Data.Response/MaxResponse),10)
    %ylim([0,max(Data.CPM)]);
    xlabel('Orientation [Deg]');ylabel('CPM');zlabel('Response [Spikes/Sec]')
    title('1st Stimulus Switch');
    ax2 = nexttile;
    contourf(Data.Orientations,Data.CPM,(Data.AltResponse/MaxResponse),10)
    %ylim([0,max(Data.CPM)]);
    xlabel('Orientation [Deg]');ylabel('CPM');zlabel('Norm Response')
    title('2nd Stimulus Switch');
    cb = colorbar;
    cb.Layout.Tile = 'east';
    clim(ax1,[0 1]); clim(ax2,[0 1]);
    % Plot Response per Strip Width
    StripWidths = 1000./(Data.CPM*2);

    CombinedResponse = (Data.Response+Data.AltResponse)/2;
    % Find the maximum value in the matrix and its linear index
    [maxValue, linearIndex] = max(CombinedResponse(:));
    % Convert the linear index to row and column indices
    [MaxCPM, MaxOR] = ind2sub(size(CombinedResponse), linearIndex);

    figure();
    plot(StripWidths,Data.Response(:,MaxOR)/Data.Response(MaxCPM, MaxOR))
    xlabel('Grating Strip Width [\mum]')
    ylabel('Normilized Firing Rate')
    xlim([0 200]); ylim([0 1.1]);
end
%% Natural Vis Intensity Curve
RateCalcWindow = binsize_sec*1000; % PSTH Bin Size in ms
if StimTime > 0.5
    SponFlag2 = 1;
else
    SponFlag2 = 0;
end
[Data] = NaturalIntensityCalc(Data,RateCalcWindow,fname,WindowSize,Psth_sort(Data.Clusters,:),sampling_freq,SponFlag2);
%% Prosthetic Intensity Response Curve
RateCalcWindow = binsize_sec*1000; % PSTH Bin Size in ms
[Data] = ProstheticIntensityCalc(Data,RateCalcWindow,fname,WindowSize,Psth_sort(Data.Clusters,:));
%% CPD Selectivity over different data files
RateCalcWindow = binsize_sec*1000; % PSTH Bin Size in ms
[Data] = CPDCalc(Data,RateCalcWindow,fname,sampling_freq,Rast_sort{Data.Clusters},Psth_sort(Data.Clusters,:));

%% Save Data variable
% Change "SavePath" to the desiered folder
%SavePath ='C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Results\Data\CFF\';
%SavePath ='C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Results\Data\MicronAlt\';
SavePath ='C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Proccesed - Data\';
%SavePath = 'C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Results\Data\MicronHDMICPD\';
stimfreq = num2str(round(1/mean(diff(stimulus_times))));
[a,b] = regexp(fname,'_\d*ms');
if ~isempty(a)
    Data.StimDur = fname(a+1:b-2);
end
[a,b] = regexp(pathname,'Data\\\S{10,10}');
date = pathname(a+5:b); date = datestr(datetime(strrep(date,'.','-')),'dd-mm-yyyy');
[a,b] = regexp(fname,'Loc\d*');Data.Loc = fname(a+3:b);
%Ch = cell2mat(inputdlg('Insert the Number of the Channel'));
Ch = answer{3};
%save([SavePath,'AltNatural_Loc',Data.Loc,'_','Ch',Ch,'_',date,'.mat'],"Data");
%save([SavePath,'AltProsthetic_Loc',Data.Loc,'_','Ch',Ch,'_',date,'.mat'],"Data");
%save([SavePath,date,'_FFNatural_Loc',Data.Loc,'_','Ch',Ch,'_',Data.StimDur,'ms','.mat'],"Data");
%save([SavePath,date,'_FFProsthetic_Loc',Data.Loc,'_','Ch',Ch,'_',Data.StimDur,'ms','.mat'],"Data");
save([SavePath,'CFFNatural_Loc',Data.Loc,'_','Ch',Ch,'_',Data.StimDur,'ms_',date,'.mat'],"Data");
%save([SavePath,'CFFProsthetic_Loc',Data.Loc,'_','Ch',Ch,'_',Data.StimDur,'ms_',date,'.mat'],"Data");
%save([SavePath,'DSGHDMINatural_Loc',Data.Loc,'_','Ch',Ch,'_',date,'.mat'],"Data");
%% Spectrogram
% FullRast = reshape(Rast,[1,size(Rast,1)*size(Rast,2)]);
% [S,F,T,P] = spectrogram(raw_data{1}(stimulus_indexes(1):end),1100,[],22000,sampling_freq);
% %[S,F,T,P] = spectrogram(full(FullRast),2200,[],22000,sampling_freq);
% figure();
% imagesc(T,F,10*log10(P));
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% ylim([0 20])

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
RastCount = 0;
for l=1:length(ind_rast)
    RastCount = RastCount+length(ind_rast{l});
end
%% Split Waveforms by latency
for c = 1:Data.dim{1}
    latency  = 50; % latency cutoff in ms.
    WavformSplit(latency,Data.Rast_sort{c},Data.ind_rast_sort{c},...
        Data.ClusteredspikeIdx{c},Data.SortedSpikes{1}{c},stimulus_indexes,sampling_freq)
    % title(['Cluster ', num2str(c),' Latency: ',num2str(latency)])

end

