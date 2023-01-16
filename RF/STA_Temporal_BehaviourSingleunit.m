function [STA_cent,STA_Surr]=STA_Temporal_BehaviourSingleunit(e_center,e_surr1,e_surr2,Frame,AC)


%% mean of pixels in this area... for all 10 frames

for j=1:size(Frame,1)
%  for j=33
    for k=1:size(Frame,2)
%  for k=2
        for m=10:-1:1
%        
          

% M=Frame{j,k}(50:end,90:end,m)-mean(mean(Frame{j,k}(50:end,90:end,m)));
M=Frame{j,k}(50:end,90:end,m)-mean(mean(Frame{j,k}(50:end,90:end,m)));

[row, col]=find(e_center{j,k});
% [row_surr, col_surr]=find(e_surr1{j,k}-e_center{j,k});
% [row_surr, col_surr]=find(e_surr2{j,k}-e_center{j,k});
[row_surr, col_surr]=find(e_surr1{j,k}-e_surr2{j,k});

% [row_surr, col_surr]=find(e_surr2{j,k});

%         
%            STA_cent{j,k}(m)=mean(mean(M(min(row):max(row),min(col):max(col))))-mean(mean(M));
           STA_cent{j,k}(m)=mean(mean(M(min(row):max(row),min(col):max(col))));

%  STA_Surr{j,k}(m)=mean(mean(M(min(row_surr):max(row_surr),min(col_surr):max(col_surr))))-mean(mean(M));
  STA_Surr{j,k}(m)=mean(mean(M(min(row_surr):max(row_surr),min(col_surr):max(col_surr))));

        end
    end
end
%% Normalize Values Between 0-1
for i = 1:length(STA_cent)
    STA_cent{i} = STA_cent{i} / max(STA_cent{i});
    STA_Surr{i} = STA_Surr{i} / max(STA_Surr{i});
end
%% plot
t_sta=-450:50:0;

%ChannelPosition = [1,9,5,6,14,13,2,10,15,7,12,11,3,4,16,8];
ChannelPosition = [1,2];
%AC = [1:16];
indx_elec=[];
%Channels_Psth= figure('name','Receptive Field');
% count=1;
% range = max(max(max(cell2mat(Frame))));
for k=1:size(Frame,2)
    
    for i=1:length(ChannelPosition)
    figure;
        %     figure('Name',['Elec:',num2str(i)]);


        count = ChannelPosition(i);

%         subplot(1,2,count)
        indx_elec =i;

%         plot(t_sta(end:-1:1),(STA_Surr{indx_elec(end),k}),'b')  
%         hold on
        plot(t_sta(end:-1:1),( STA_cent{indx_elec(end),k}),'r')
        yticklabels
        title(['Channel: ',num2str(AC(i))]);
        ylim([min(STA_cent{indx_elec(end),k})-0.1 max(STA_cent{indx_elec(end),k})+0.1]);

    end
end









 
 
 
 
% plot(t(end:-1:1),st_cent-st_cent(10),'r')
% hold on
% plot(t(end:-1:1),st_surr-st_surr(10),'b')