function [aligned_waveforms, waveforms]=spike_waves29_11(spike_ind,fs,chan_filt,std_factor)

for i=1:length(spike_ind)
    for m=1:length(spike_ind{i})
        if spike_ind{i}(m)-5*10^-3*fs>0 &spike_ind{i}(m)+5*10^-3*fs<length(chan_filt{i});
        waveforms{i}(:,m)=chan_filt{i}(spike_ind{i}(m)-5*10^-3*fs:spike_ind{i}(m)+5*10^-3*fs);    
        end
    end
end
%% align spikes 

for i=1:length(spike_ind)
 if ~isempty(waveforms{i})
    for m=1:size(waveforms{i},2)
        if std_factor>0
            [val ind]=max(waveforms{i}(:,m));
        else 
                        [val ind]=min(waveforms{i}(:,m));
        end
          if ind-1.5*10^-3*fs>0 &ind+1.5*10^-3*fs<size(waveforms{i},1);
        aligned_waveforms{i}(:,m)=waveforms{i}(ind-1.5*10^-3*fs:ind+1.5*10^-3*fs,m);
    
          end
    end
end
end

            