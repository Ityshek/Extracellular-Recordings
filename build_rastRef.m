function [Rast,Spike,Av_spike,events_ind,ind_rast,spike_stim,spike_Times]=build_rastRef(stimulus_indexes,Stim_times,raw_data,sampling_freq,thresh,outlier)


Stim_freq=round((1/mode(diff(Stim_times))),2);
%Stim_freq = 1; % for 1Hz 10 frame trigger
ISI=round(median(diff(Stim_times(2:end))),2);
%ISI=1; % for 1Hz 10 frame trigger
Rast=sparse(round(length(stimulus_indexes)/ceil(Stim_freq)),round(ISI*sampling_freq+10*10^-3*sampling_freq)); % for 1Hz trigger
%Rast=sparse(round(length(stimulus_indexes)/ceil(Stim_freq)),round(mean(diff(Stim_times)),2)*sampling_freq+10*10^-3*sampling_freq); % for 1Hz 10 frame trigger

count_stim=1;



events_ind=find(circshift(raw_data,[0 1])>thresh&raw_data<thresh); % indecies for spike times, detected by thresh.
% while min(diff(events_ind))<0.003*sampling_freq
%     for i=2:length(events_ind)
%         if ~isempty(events_ind(i-1))
%             if events_ind(i) < events_ind(i-1)+0.003*sampling_freq
%                 events_ind(i) = NaN;
%             end
%         end
%     end
%     events_ind = rmmissing(events_ind);
% end
for i=1:length(events_ind)
    if events_ind(i)-1*10^-3*sampling_freq>0&events_ind(i)+1.5*10^-3*sampling_freq<length(raw_data) %&circshift(raw_data(events_ind(i)),[0 -1])>-200
        if max([raw_data(events_ind(i)-((1*10^-3)*sampling_freq):events_ind(i)+((1.5*10^-3)*sampling_freq))]) < outlier
            spike{i}=raw_data(events_ind(i)-1*10^-3*sampling_freq:events_ind(i)+1.5*10^-3*sampling_freq); % saves the  i spike's amplitudes
        end
    end
end


%% build raster
%for k=1:2:length(stimulus_indexes) % loop over 1 every 2 triggers.
for k=1:length(stimulus_indexes)
    counter_events=1;
    for i=1:length(events_ind) %loop over all events(spikes found)



        if events_ind(i)>(stimulus_indexes(k)-(10*10^-3)*sampling_freq) && events_ind(i)<stimulus_indexes(k)+(ISI-10^-3)*sampling_freq
            if max([raw_data(events_ind(i)-((1*10^-3)*sampling_freq):events_ind(i)+((1.5*10^-3)*sampling_freq))]) < outlier

                ind_rast{count_stim}(counter_events)=events_ind(i)-(stimulus_indexes(k)-(10*10^-3)*sampling_freq);
                spike_stim{count_stim}(:,counter_events)=raw_data(events_ind(i)-1*10^-3*sampling_freq:events_ind(i)+1.5*10^-3*sampling_freq); %contains 3 msec around neg peak of each spike.
                counter_events=counter_events+1;
            else
                events_ind(i) = NaN;
                spike{i} = NaN;
            end

        end

    end
    count_stim=count_stim+1;
end

events_ind = rmmissing(events_ind);
count = 1;
for k=1:length(spike)
    if ~isnan(spike{k})
        Spike{count} = spike{k};
        count = count+1;
    end
end

EmptyCheck = exist('ind_rast','var');
if EmptyCheck ~=0
    for m=1:length(ind_rast)
        if ~isempty((ind_rast{m}))

            Rast(m,round(ind_rast{m}))=1;


        end
    end
    Av_spike=mean(cell2mat(spike),1);
    spike_Times = events_ind/sampling_freq;
else
    ind_rast = [];
    spike_stim{count_stim} = [];
    Av_spike = [];
    spike = [];
    spike_Times = [];
    indx_spike = [];
end
%%



