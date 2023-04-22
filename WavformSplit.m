function [Data] = WavformSplit(latency,Rast_sort,ind_rast_sort,idxspike,spike,stimulus_indexes,sampling_freq)
idxlatency = latency/1000*sampling_freq;
ArtifactGroup = []; SpikeGroup = [];
idxspike = rmmissing(idxspike);
for i = 1:length(ind_rast_sort)
    [spikes,ids] = intersect(idxspike,ind_rast_sort{i}+stimulus_indexes(i)-440);
    if length(spikes) > 0
        for k = 1:length(ids)
            if idxspike(ids(k)) < stimulus_indexes(i) + idxlatency
                ArtifactGroup = [ArtifactGroup;spike(:,ids(k))']; % Suspected Artifact
            else
                SpikeGroup = [SpikeGroup;spike(:,ids(k))']; % Suspected Clean
            end
        end
    end
end
%% Plot Waveforms of the groups
spikedur = round(size(SpikeGroup,2)/sampling_freq*1000,1); % duration of spike in ms
CenterIdx = floor((spikedur/2)*(sampling_freq/1000));
sizems = 2;
WaveSize = sizems*(sampling_freq/1000); 
t_spike = linspace(0,sizems,WaveSize+1);

Groups{1} = ArtifactGroup; Groups{2} = SpikeGroup;
Titles = {'Suspected Artifact','Suspected Clean'};
figure();
for i = 1:length(Groups)
    for k = 1:size(Groups{i},1)
        [CenteredSpikeVal,CenteredSpikeIdx] = min(Groups{i}(k,:));
        if CenteredSpikeIdx < 178 && CenteredSpikeIdx > 44
        CenteredSpike(k,:) = Groups{i}(k,CenteredSpikeIdx-WaveSize/2:CenteredSpikeIdx+WaveSize/2);
        subplot(1,2,i)
        plot(t_spike,CenteredSpike(k,:))
        hold on
        end
    end
    plot(t_spike,mean(CenteredSpike,1),"LineWidth",2,"Color","k")
    xlabel('Time [ms]')
    ylabel('Amplitude [\muV]')
    title(Titles{i})
    ylim([-150 150])
end
%% Feature Comparison

P2P = cell(1,length(Groups)); T2P = cell(1,length(Groups)); T2Z = cell(1,length(Groups));
for g = 1:length(Groups)
    for i = 1:size(Groups{g},1)
        [minval,minidx] = min(Groups{g}(i,:));
        [maxval,maxidx] = max(Groups{g}(i,minidx:end));
        maxidx = maxidx + minidx-1;
        t2z = min(find(Groups{g}(i,maxidx:end) ==0)/sampling_freq);
        if ~isempty(t2z)
            T2Z{g} = [T2Z{g};t2z];
            P2P{g} = [P2P{g};abs(maxval)/abs(minval)];
            T2P{g} = [T2P{g};(maxidx-minidx)/sampling_freq];
        end
    end
end

% Plot 
figure();
col = ['k','b'];
for g = 1:length(Groups)
plot3(P2P{g},T2P{g},T2Z{g},'*',Color=col(g))
hold on
end
ylabel('Time to Max Peak [ms]')
xlabel('Max to Min Ratio')
zlabel('Time to Zero [ms]')
legend('Artifacts', 'Spikes')
