function [Rast_sort]=build_rast_sort2(idx,stimulus_indexes,sampling_freq,ind_rast_re2,dim,Stim_times)
% [Amplitude, Width,Spikes]=Sort_spikes(raw_data,sampling_freq);
% [width1,width2,width3,idx,Spikes1,Spikes2,Spikes3]=cluster_data(Amplitudes, Width,Spikes);
ISI=round(mean(diff(Stim_times(2:end))),2);
for i=1:dim
    Rast_sort{i}=sparse(length(ind_rast_re2),ISI*sampling_freq+10*10^-3*sampling_freq);
end

for i=1:length(ind_rast_re2)    
        k=length(ind_rast_re2{i});
        idx_2{i}=idx(i:i+k-1);
end

for k=1:dim
    for i=1:length(idx_2)
        Rast_sort{k}(i,ind_rast_re2{i}(idx_2{i}==k))=1;
    end
end







% thresh=std_Factor*std(double(raw_data));
% count=1;
% for k=1:length(stimulus_indexes)
%     if stimulus_indexes(k)-1*10^-3*sampling_freq>0&stimulus_indexes(k)+1000*10^-3*sampling_freq<length(raw_data)
%         %
%         temp_seg=raw_data(stimulus_indexes(k)-1*10^-3*sampling_freq:stimulus_indexes(k)+1000*10^-3*sampling_freq);
%             ind_rast{k}=find(circshift(temp_seg,[0 1])>thresh&temp_seg<thresh);
% %         [spike_indexes1] =find_Spikes(temp_seg,sampling_freq);
% %         ind_rast{k}=spike_indexes1;
%         ind_rast_re{k}=ind_rast{k}(find(diff(ind_rast{k})>sampling_freq*0.002));
%         for i=1:length(ind_rast_re{k})
%             Spikes{k}{i}=raw_data(ind_rast_re{k}(i)+stimulus_indexes(k)-1.5*10^-3*sampling_freq:ind_rast_re{k}(i)+stimulus_indexes(k)+1.5*10^-3*sampling_freq);
%             [val_max(i), ind_max(i)]= max(Spikes{k}{i});
%             [val_min(i), ind_min(i)]=min(Spikes{k}{i});
%             Amplitude{k}(i)=val_max(i)-val_min(i);
%             Width{k}(i)=abs(ind_max(i)-ind_min(i))/sampling_freq;
%             indx_spike(count)=ind_rast_re{k}(i)+round(stimulus_indexes(k));
%             count=count+1;
%         end
%
%     end
% end
%
%
% %% sort the indexes
% for k=1:length(stimulus_indexes)
%     count1=1;
%     count2=1;
%     count3=1;
%     [amplitude1{k},amplitude2{k},amplitude3{k},idx{k},Spikes1{k},Spikes2{k},Spikes3{k}]=cluster_data(Amplitude{k}, Width{k}, Spikes{k});
%     for m=1:length(idx{k})
%         if idx{k}(m)==1
%             ind_rast1{k}(count1)=ind_rast_re{k}(m);
%             Spikes_sort1{k}{count1}=Spikes{k}{m};
%             width1{k}(count1)=Width{k}(m);
%             count1=count1+1;
%
%         elseif  idx{k}(m)==2
%             ind_rast2{k}(count2)=ind_rast_re{k}(m);
%             Spikes_sort2{k}{count2}=Spikes{k}{m};
%             width2{k}(count2)=Width{k}(m);
%
%             count2=count2+1;
%
%
%         elseif idx{k}(m)==3
%             ind_rast3{k}(count3)=ind_rast_re{k}(m);
%             Spikes_sort3{k}{count3}=Spikes{k}{m};
%             width3{k}(count3)=Width{k}(m);
%             count3=count3+1;
%
%
%         end
%
%     end
%
% end
%
% %% build rast
% if ~isempty(ind_rast1)
%     for k=1:length(ind_rast1)
%      Rast_sort{1}(k,ind_rast1{k})=1;
%     end
% else
%     Rast_sort{1}=[];
%     width1=[];
% end
%
% if ~isempty(ind_rast2)
%     for k=1:length(ind_rast2)
%      Rast_sort{2}(k,ind_rast2{k})=1;
%     end
% else
%     Rast_sort{2}=[];
%     width2=[];
% end
%
% if ~isempty(ind_rast3)
%     for k=1:length(ind_rast3)
%      Rast_sort{3}(k,ind_rast3{k})=1;
%     end
% else
%     Rast_sort{2}=[];
%     width3=[];
% end
%
%
