function Sorted_Spikes=sort_spikes3(idx,Aligned_Spikes,dim)
Sorted_Spikes = {};
%
% for k=1:3
%     for i=1:length(idx)
%         Sorted_Spikes{k}=nan(size(Aligned_Spikes{i}));
%     end
% end
if ~isempty(idx)
for k=1:dim
    count=1;
    for i=1:length(idx)
        %         for m=1:size(Aligned_Spikes{i},2)
        if idx(i)==k
            Sorted_Spikes{k}(:,count)=Aligned_Spikes(i,:);
            count=count+1;
        end
    end
end
end