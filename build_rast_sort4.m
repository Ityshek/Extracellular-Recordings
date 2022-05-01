function [Rast_sort]=build_rast_sort4(idx,stimulus_indexes,sampling_freq,dim,Aligned_idx,Stim_times)
% [Amplitude, Width,Spikes]=Sort_spikes(raw_data,sampling_freq);
% [width1,width2,width3,idx,Spikes1,Spikes2,Spikes3]=cluster_data(Amplitudes, Width,Spikes);
StimDur = int32(mean(diff(stimulus_indexes)));
for i=1:dim
    Rast_sort{i}=sparse(length(stimulus_indexes),StimDur+1*10^-3*sampling_freq);
end
count_stim=1;
ISI=round(mean(diff(Stim_times(2:end))),2);
counter_events2=1;
for k=1:length(stimulus_indexes)
    counter_events=1;
    for i = 1:length(Aligned_idx)
        if Aligned_idx(i)>stimulus_indexes(k)-10*10^-3*sampling_freq && Aligned_idx(i)<stimulus_indexes(k)+(ISI-10^-3)*sampling_freq
            if Aligned_idx(i)-stimulus_indexes(k)>0
                ind_rast_re2{count_stim}(counter_events)=Aligned_idx(i)-stimulus_indexes(k);
                UpdateIdx(counter_events2) = idx(i);
                counter_events=counter_events+1;
                counter_events2=counter_events2+1;
            end
        end
    end
    count_stim=count_stim+1;
end

count = 1;
for i=1:length(ind_rast_re2)
    k=length(ind_rast_re2{i});
    idx_2{i}=UpdateIdx(count:count+k-1);
    count = count+k;
end

for k=1:dim
    for i=1:length(idx_2)
        Rast_sort{k}(i,ind_rast_re2{i}(idx_2{i}==k))=1;
    end
end







