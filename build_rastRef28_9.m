function [Rast,spike,Av_spike,events_ind,ind_rast,spike_stim,spike_Times]=build_rastRef28_9(stimulus_indexes,raw_data,sampling_freq,thresh,IdxISI)


Stim_times=stimulus_indexes/sampling_freq;
Stim_freq=round((1/mean(diff(Stim_times(1:end)))),2);
ISI=round(mean(diff(Stim_times(2:end))),2);
Rast=sparse(round(length(stimulus_indexes)/ceil(Stim_freq)),int64(ISI*sampling_freq)+10*10^-3*sampling_freq);
count_stim=1;
spike = []; Av_spike = []; spike_stim = []; spike_Times = [];

events_ind=find(circshift(raw_data,[0 -1])>thresh&raw_data<thresh); % indecies for spike times, detected by thresh.
% while min(diff(events_ind))<0.002*sampling_freq
%     for i=2:length(events_ind)
%         if ~isempty(events_ind(i-1))
%             if events_ind(i) < events_ind(i-1)+0.002*sampling_freq
%                 events_ind(i) = NaN;
%             end
%         end
%     end
% events_ind = rmmissing(events_ind);    
% end
for i=1:length(events_ind)
    if events_ind(i)-0.5*10^-3*sampling_freq>0&events_ind(i)+1.5*10^-3*sampling_freq<length(raw_data)%&circshift(raw_data(events_ind(i)),[0 -1])>-200
        spike{i}=raw_data(events_ind(i)-1*10^-3*sampling_freq:events_ind(i)+1*10^-3*sampling_freq); % saves the  i spike's amplitudes
    end
end


%% build raster

for k=1:length(stimulus_indexes) % loop over all triggers from the 2nd to last.
    counter_events=1;
    ind_rast{count_stim} = [];
    for i=1:length(events_ind) %loop over all events(spikes found)
        if events_ind(i)>stimulus_indexes(k)-10*10^-3*sampling_freq && events_ind(i)<stimulus_indexes(k)+(ISI-10^-3)*sampling_freq
            if events_ind(i)-stimulus_indexes(k)>0
                ind_rast{count_stim}(counter_events)=events_ind(i)-stimulus_indexes(k);
                spike_stim{count_stim}(:,counter_events)=raw_data(events_ind(i)-0.5*10^-3*sampling_freq:events_ind(i)+1.5*10^-3*sampling_freq); %contains 1 sec in each cell, and all spikes in that 1 sec within each cell.
                counter_events=counter_events+1;
%             elseif exist("ind_rast","var")
%                 ind_rast{count_stim}(counter_events) = [];
            end           
        end
    end
    count_stim=count_stim+1;
end



for m=1:length(ind_rast)
    if ~isempty(ind_rast{m})
        
        Rast(m,ind_rast{m})=1;
        
    end
end
%%
if ~isempty(spike)
Av_spike=mean(cell2mat(spike),1);
end
spike_Times = events_ind/sampling_freq;


