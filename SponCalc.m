function [Spon,SponStd] = SponCalc(Data,sampling_freq,stimulus_indexes,window)

SponVec = zeros(1,length(stimulus_indexes)-1);
Window = (window(2)-window(1))*0.001; % Define size of window to calc spikes in

for k = 1:length(Data.Clusters)
    for i = 2:length(stimulus_indexes)
        A = find((Data.ClusteredspikeIdx{Data.Clusters(k)} > stimulus_indexes(i)-sampling_freq*Window));
        B =  find(Data.ClusteredspikeIdx{Data.Clusters(k)} < stimulus_indexes(i));
        SponVec(i-1) = length(intersect(A,B));
    end
    %Data.Spon{k} = [Data.Spon{k};mean(SponVec)];
    %Data.SponStd{k} = [Data.SponStd{k};std(SponVec)];

    % Calc Spon for with smoothed avarege
if stimulus_indexes(1) >= 30*sampling_freq % Check if there are at least 30s of idle recording 
NumSpikes = length(find(Data.ClusteredspikeIdx{Data.Clusters(k)} < stimulus_indexes(1))); % Look for all spikes prior to 1st trigger
LatSponVec = zeros(1,stimulus_indexes(1)); % Build spike train vector
B = find(~isnan(Data.ClusteredspikeIdx{Data.Clusters(k)}(1:NumSpikes)));
LatSponVec(Data.ClusteredspikeIdx{Data.Clusters(k)}(B)) = 1; % insert spikes to vector
LatencyWindow = 0.05*sampling_freq; % Define window to calc spikes, should be same as smoothing of PSTH.
%Win = ones(1,LatencyWindow);
SmoothedSpon = smooth(LatSponVec,LatencyWindow);
Spon{k} = mean(SmoothedSpon)*window(1)*10^-3*sampling_freq; % Calc mean spon in same bin size as PSTH.
SponStd{k} = std(SmoothedSpon)*window(1)*10^-3*sampling_freq;
else % calc spon over second half of the PSTH
    Spon{k} = mean(Data.Psth_sort{end}(end/2+1:end));
    SponStd{k} =std(Data.Psth_sort{end}(end/2+1:end));
end
end