function [Data] = SponCalc(Data,sampling_freq,stimulus_indexes,window)

SponVec = zeros(1,length(stimulus_indexes)-1);
Window = (window(2)-window(1))*0.001;

for k = 1:length(Data.Clusters)
    for i = 2:length(stimulus_indexes)
        A = find((Data.ClusteredspikeIdx{Data.Clusters(k)} > stimulus_indexes(i)-sampling_freq*Window));
        B =  find(Data.ClusteredspikeIdx{Data.Clusters(k)} < stimulus_indexes(i));
        SponVec(i-1) = length(intersect(A,B));
    end
    Data.Spon{k} = [Data.Spon{k};mean(SponVec)];
    Data.SponStd{k} = [Data.SponStd{k};std(SponVec)];
end
end