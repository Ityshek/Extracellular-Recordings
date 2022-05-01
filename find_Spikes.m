function [spike_indexes1 ,spike_indexes] =find_Spikes(raw_data,sampling_freq,std_Factor)
Event_threshold=std_Factor*std(raw_data);
% spike_indexes=find((circshift(raw_data,[0 1])>Event_threshold)&&raw_data<Event_threshold);
% spike_indexes1=find((circshift(raw_data,[0 1])<Event_threshold)&raw_data>Event_threshold);
% 
spike_indexes1=find(raw_data>Event_threshold);
temp_indx=find((diff(spike_indexes1)/sampling_freq)<0.003);
spike_indexes1(temp_indx)=0;
spike_indexes1=spike_indexes1(find((spike_indexes1)));
for i=1:length(spike_indexes1)
     [~,max_val(i)]=max(raw_data(spike_indexes1(i)-1.5*10^-3*sampling_freq:spike_indexes1(i)+1.5*10^-3*sampling_freq));
 spike_indexes(i)=spike_indexes1(i)+max_val(i)-1.5*10^-3*sampling_freq;
end
%  