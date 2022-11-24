function [Spikes_Matrix, Average_Spike,Aligned_idx]=Align_spikes4(spike_stim,fs,std_Factor,indx_spike)
Aligned_idx = nan(1,length(spike_stim));
for k=1:length(spike_stim)
    if isempty(spike_stim{1,k})~=1
        
        if std_Factor>0
            [max_val(k),max_index(k)]=max(spike_stim{1,k});
        else
            [max_val(k),max_index(k)]=min(spike_stim{1,k});
        end
        spike_indexes(k)=max_index(k);
        if spike_indexes(k)+2*10^-3*fs<=length(spike_stim{1,k})&&spike_indexes(k)-1.5*10^-3*fs>0   
            Aligned_Spikes{k}=spike_stim{1,k}(spike_indexes(k)-1.5*10^-3*fs:spike_indexes(k)+2*10^-3*fs);
            Aligned_idx(k) = indx_spike(k);
        end
    end
end
count = 1;
for i=1:length(Aligned_Spikes)
    if size(cell2mat(Aligned_Spikes(i)),2)>1
        Spikes_Matrix(count,:) = cell2mat(Aligned_Spikes(i));
        count = count+1;
    end
end
Average_Spike=mean(Spikes_Matrix);
Aligned_idx = rmmissing(Aligned_idx);