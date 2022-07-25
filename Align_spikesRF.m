function [Spikes_Matrix, Average_Spike,Aligned_idx]=Align_spikesRF(AC,spike_stim,fs,std_Factor,indx_spike)

for c=1:length(AC)
    Aligned_Spikes = [];
    Aligned_idx{AC(c)} = nan(1,length(spike_stim{AC(c)}));
for k=1:length(spike_stim{AC(c)})
    if isempty(spike_stim{AC(c)}{1,k})~=1
        
        if std_Factor>0
            [max_val(k),max_index(k)]=max(spike_stim{AC(c)}{1,k});
        else
            [max_val(k),max_index(k)]=min(spike_stim{AC(c)}{1,k});
        end
        spike_indexes(k)=max_index(k);
        if spike_indexes(k)+1.5*10^-3*fs<=length(spike_stim{AC(c)}{1,k})&&spike_indexes(k)-1*10^-3*fs>0   
            Aligned_Spikes{k}=spike_stim{AC(c)}{1,k}(spike_indexes(k)-1*10^-3*fs:spike_indexes(k)+1.5*10^-3*fs);
            Aligned_idx{AC(c)}(k) = indx_spike{AC(c)}(k);
        end
    end
end
count = 1;
if ~isempty(Aligned_Spikes)
for i=1:length(Aligned_Spikes)
    if size(cell2mat(Aligned_Spikes(i)),2)>1
        Spikes_Matrix{AC(c)}(count,:) = cell2mat(Aligned_Spikes(i));
        count = count+1;
    end
end
Average_Spike{AC(c)}=mean(Spikes_Matrix{AC(c)});
Aligned_idx{AC(c)} = rmmissing(Aligned_idx{AC(c)});
aligned_spikes{AC(c)} = Aligned_Spikes;
else
Average_Spike{AC(c)}=[];
Aligned_idx{AC(c)} = [];
aligned_spikes{AC(c)} = [];
end    
end