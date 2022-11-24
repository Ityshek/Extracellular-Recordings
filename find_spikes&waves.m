function [Events_ind,spike]=find_spikes&waves(AC,raw_data,sampling_freq,std_factor,stim_indx)
% global thresh
% global std_factor
outlier = 400;
Count = 1;
for c=1:length(AC)
%     std_factor=10;
    count=1;
    % detect spontaneous  events based on threshoding using stdv
    %     thresh(i)=std_factor*std(chan_filt{i});%thresholding ...
    thresh(AC(i))=std_factor*nanstd(clean_chan{AC(i)});%thresholding ...
%     if thresh(i)>20
    temp_seg=clean_chan{AC(i)}(1:stim_indx{AC(i)}(end)+fs*0.1); % takes the raw data up to 100ms after the last frame presented
    spike_ind{AC(i)}=find(circshift(temp_seg,[0 1])>thresh(AC(i))&temp_seg<thresh(AC(i)));
    %             ind_rast_re{i}{count}=ind_rast{i}{k}(find(diff(ind_rast{i}{k})>fs*0.01));
    spike_ind{AC(i)}=spike_ind{AC(i)}(find(diff(spike_ind{AC(i)})>=(0.003*fs)));
%     else
%         spike_ind{i}=[];
%     end



Events_ind=find(circshift(raw_data{AC(c)},[0 1])>thresh&raw_data{AC(c)}<thresh); % indecies for spike times, detected by thresh.

% remove events inside refractory period
events_ind{AC(c)} = [];
for i=2:length(Events_ind)
    if Events_ind(i)+(1.5*10^-3)*sampling_freq<=length(raw_data{AC(c)}) & Events_ind(i)-(1.5*10^-3)*sampling_freq>0 
    if Events_ind(i) - Events_ind(i-1) > 88 & (max([raw_data{AC(c)}(Events_ind(i)-((1.5*10^-3)*sampling_freq):Events_ind(i)+((1.5*10^-3)*sampling_freq))]) < outlier)
        events_ind{AC(c)} = [events_ind{AC(c)},Events_ind(i)];
    end
end
end

for i=1:length(events_ind{AC(c)})
    if events_ind{AC(c)}(i)-4*10^-3*sampling_freq>0&events_ind{AC(c)}(i)+4*10^-3*sampling_freq<length(raw_data{AC(c)}) %&circshift(raw_data(events_ind(i)),[0 -1])>-200
        if max([raw_data{AC(c)}(events_ind{AC(c)}(i)-((1*10^-3)*sampling_freq):events_ind{AC(c)}(i)+((1.5*10^-3)*sampling_freq))]) < outlier
            spike{AC(c)}{i}=raw_data{AC(c)}(events_ind{AC(c)}(i)-4*10^-3*sampling_freq:events_ind{AC(c)}(i)+4*10^-3*sampling_freq); % saves the  i spike's amplitudes
        end
    end
end
Count = Count+1;
end