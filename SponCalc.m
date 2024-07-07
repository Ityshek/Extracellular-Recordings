function [Spon,SponStd] = SponCalc(PsthBinned,Rast_sort,sampling_freq,stimulus_indexes,window,raw_data,SponRelevantSpikes,Flag,RateWinSize,WindowSize)
if Flag == '0'
    %% Calc From Period Between Stimuli

    %BetweenWindow = (window(2)-window(1))*0.01*sampling_freq;
    %SponMat = full(Rast_sort(:,end-BetweenWindow-2200:end-2201)); % Extract the baseline from the rasters in a window prior to each stimulus.
    %     BetweenWindow = RateWinSize*sampling_freq;
    %     SponMat = full(Rast_sort(:,end-BetweenWindow-2200:end-2201)); % Extract the baseline from the rasters in a window prior to each stimulus.
    %     Spon = sum(sum(SponMat))/size(SponMat,1)*(sampling_freq/size(SponMat,2)); % Calc average Baseline in spikes/sec
    %     SponStd = std(sum(SponMat,2))*(sampling_freq/size(SponMat,2)); % Calc std of Baseline in spikes/sec
    %
    BetweenWindow = (WindowSize/1000)/RateWinSize;
    if iscell(PsthBinned)
        for p = 1:length(PsthBinned)
            SponMat(p,:) = PsthBinned{p}(:,end-BetweenWindow:end);
        end
        Spon = mean2(SponMat); SponStd = std2(SponMat);
    else
       Spon = mean(PsthBinned(:,end-BetweenWindow:end));
       SponStd = std(PsthBinned(:,end-BetweenWindow:end));
    end
else
    %% Calc from Period Prior to 1st Stimulus
    %BinSize = (window(2) - window(1))*0.01*sampling_freq; % define binsize to clac spikes in
    BinSize = RateWinSize*sampling_freq;
    count=1;
    SponRelevantSpikes = SponRelevantSpikes(find(SponRelevantSpikes<stimulus_indexes(1))); % Pick only spikes prior to 1st stim
    for i=1:BinSize:stimulus_indexes(1)-BinSize
        SponVec(count) = sum(i<SponRelevantSpikes & SponRelevantSpikes<i+BinSize);
        count=count+1;
    end
    RateMultiplier = sampling_freq/BinSize;
    Spon = mean(SponVec*RateMultiplier); % Calc average baseline spiking rate in spikes/sec
    SponStd = std(SponVec*RateMultiplier); % Calc std of baseline spiking rate in spikes/sec
end
%%
% SponVec = zeros(1,length(stimulus_indexes)-1);
% Window = (window(2)-window(1))*0.001; % Define size of window to calc spikes in
%
% for k = 1:length(Data.Clusters)
%     for i = 2:length(stimulus_indexes)
%         A = find((Data.ClusteredspikeIdx{Data.Clusters(k)} > stimulus_indexes(i)-sampling_freq*Window));
%         B =  find(Data.ClusteredspikeIdx{Data.Clusters(k)} < stimulus_indexes(i));
%         SponVec(i-1) = length(intersect(A,B));
%     end
%     %Data.Spon{k} = [Data.Spon{k};mean(SponVec)];
%     %Data.SponStd{k} = [Data.SponStd{k};std(SponVec)];
%
% % Calc Spon for with smoothed avarege
% if stimulus_indexes(1) >= 30*sampling_freq % Check if there are at least 30s of idle recording
% NumSpikes = length(find(Data.ClusteredspikeIdx{Data.Clusters(k)} < stimulus_indexes(1))); % Look for all spikes prior to 1st trigger
% LatSponVec = zeros(1,stimulus_indexes(1)); % Build spike train vector
% B = find(~isnan(Data.ClusteredspikeIdx{Data.Clusters(k)}(1:NumSpikes)));
% LatSponVec(Data.ClusteredspikeIdx{Data.Clusters(k)}(B)) = 1; % insert spikes to vector
% LatencyWindow = 0.01; % Define window to calc spikes, should be same as of PSTH.
% %Win = ones(1,LatencyWindow);
% %SmoothedSpon = smooth(LatSponVec,LatencyWindow);
% % Calc Noise
% Sections = 25;
% for n=1:Sections
%     NoiseMean(n) = length(find(Data.ClusteredspikeIdx{Data.Clusters(k)}<((stimulus_indexes(1)/Sections)*n)));
%     if n>1
%         NoiseMean(n) = NoiseMean(n)-sum(NoiseMean(1:n-1));
%     end
% end
% NoiseMean = NoiseMean/(stimulus_indexes(1)/Sections/sampling_freq)*LatencyWindow;
% Spon{k} = mean(NoiseMean,'omitnan'); %*window(1)*10^-3*sampling_freq; % Calc mean spon in same bin size as PSTH.
% SponStd{k} = std(NoiseMean,'omitnan'); %*window(1)*10^-3*sampling_freq;
% else % calc spon over second half of the PSTH
%     Spon{k} = mean(Data.Psth_sort{end}(end/2+1:end));
%     SponStd{k} =std(Data.Psth_sort{end}(end/2+1:end));
% end
% end