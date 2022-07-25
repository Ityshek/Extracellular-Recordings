function [events_ind,spike,thresh]=find_spikesANDwaves(AC,raw_data,sampling_freq,std_factor,stim_indx)

outlier = 400;

for c=1:length(AC)
thresh(AC(c))=std_factor*nanstd(raw_data{AC(c)});%thresholding ...

 


Events_ind=find(circshift(raw_data{AC(c)},[0 1])>thresh(AC(c))&raw_data{AC(c)}<thresh(AC(c))); % indecies for spike times, detected by thresh.

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
end
end