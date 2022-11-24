function [spike_ind]=find_spikes(AC,clean_chan,fs,std_factor,stim_indx)
% global thresh
% global std_factor
Count = 1;
for i=1:length(AC)
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
Count = Count+1;
end
