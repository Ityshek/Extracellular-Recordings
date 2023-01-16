%% load the data and find the spikes
clear all
close all
global indx stim_indx
global std_factor;
std_factor=-4.5;
ChannelPosition = [1,9,5,6,14,13,2,10,15,7,12,11,3,4,16,8];
AC = [str2num(cell2mat(inputdlg('Insert the Numbers of the Recorded Channels')))];
[raw_data, fs,stim_Data,stim_sampling_rate,Begin_record,stimulus_times,Stim_indx,FrameDuration] =load_data_ConcateMultiUnit(AC);
%clean_chan=clean_artifact(chan_filt,fs);
[Spike_ind,waveforms,thresh]=find_spikesANDwaves(AC,raw_data,fs,std_factor,Stim_indx);
[aligned_Spikes,Average_Spike,Aligned_idx]=Align_spikesRF(AC,waveforms,fs,std_factor,Spike_ind);
FrameDuration = 50;

% Raw Data
figure();
for c = 1:length(AC)
    t=[0:length(raw_data{AC(c)})-1]/fs;
    threshold = thresh(AC(c))*ones(1,length(t));
    x = nan(1,length(t));
    x(Spike_ind{AC(c)})=raw_data{AC(c)}(Spike_ind{AC(c)});
    subplot(2,length(AC)/2,c)
    plot(t,raw_data{AC(c)},'b',t,threshold,'r',t,x,'*k')
    % figure settings
    xlabel('Time[Sec]','FontSize',15)
    ylabel('Amplitude[\muV]','FontSize',15)
    title(['Channel ',num2str(c+16)])
    ylim([-200 200])
    xlim([0 max(t)])
    hold on
end
%%
% Select active channels for further analysis
AC = [str2num(cell2mat(inputdlg('Select active channels for further analysis:')))];
chan_filt = raw_data;
Aligned_Spikes = aligned_Spikes;
spike_ind = Spike_ind;
stim_indx = Stim_indx;

%% %sort data
for j=1:length(AC)
    SortingFig{j} = figure();
    if ~isempty(Aligned_Spikes{AC(j)})
        ClustEvalGAP = evalclusters(Aligned_Spikes{AC(j)},'kmeans','gap','KList',[1:3]);
        dim{AC(j)}=ClustEvalGAP.OptimalK;
%         dim{AC(j)} = 2;
        if AC(j)==1
            [ filename2, pathname_clust]=uiputfile('*.m','Select a Folder to Save Clustering Results');
        end
        indx_elec =AC(j);
        %filename_clust=[filename2(1:end-2),'_elec',num2str(indx_elec)];
        if  ~isempty(Aligned_Spikes{AC(j)})
            [idx{AC(j)},c,score]=PCA_AnalysisRF((Aligned_Spikes{AC(j)}),dim{AC(j)},SortingFig{j});
            title(['Channel: ' num2str(AC(j))])
            Sorted_temp{AC(j)}=(sort_spikes3(idx{AC(j)},(Aligned_Spikes{AC(j)}),dim{AC(j)}));
            for z=1:length(idx{AC(j)})
                for l=1:dim{AC(j)}
                    if idx{AC(j)}(z)==l

                        spike_ind_t{AC(j),l}(z)=spike_ind{AC(j)}(z);
                    end
                end
            end
        end
    end
    figure('Name',['Spike Waveforms Channel: ' num2str(AC(j))]);
    for i=1:dim{AC(j)}
        t_sort=((0:size(Sorted_temp{1,AC(j)}{i},1)-1)/fs)*10^3;
        Avg_Sorted_Spikes{AC(j)}(:,i) = mean(Sorted_temp{AC(j)}{i},2);

        subplot(1,dim{AC(j)},i)
        for u = 1:size(Sorted_temp{AC(j)}{i},2)
            plot(t_sort,Sorted_temp{AC(j)}{i}(:,u))
            hold on
        end
        plot(t_sort,Avg_Sorted_Spikes{AC(j)}(:,i),'k','linewidth',2)
        title(['Sorted waveforms - Cluster ',num2str(i)]);
        hold on
        xlabel('Time[mSec]','FontSize',20)
        ylabel('Amplitude[\muV]','FontSize',20)
        ylim([-150 100]);
    end
end

%%
[Im, ImageFiles]=load_image_directory;
for i=1:length(Im)
    %    for i=1:863

    Im_new{i}=Im{i};

    im_size1(i)=size(Im_new{i},1);
    im_size2(i)=size(Im_new{i},2);

end
%% Allocate an image for each stimulus

for i=1:length(Im_new)
    Im_new{i}=Im{i}(1:min(im_size1),1:min(im_size2));

end


for i=1:length(stim_indx{AC(1)})
    if rem(i,length(Im_new))
        im_indx{i}=Im_new{rem(i,length(Im_new))};
    else
        im_indx{i}=Im_new{length(Im_new)};

    end
end
%%
for i=1:length(AC)

    for k=1:dim{AC(i)}

        for m=1:10
            %  Frame{i}=zeros(581,601,3);
            Frame{AC(i),k}(:,:,m)=zeros(min(im_size1),min(im_size2));
        end
    end
end

%%
frame_time=stim_indx{AC(1)}/(fs);
% if isempty (frame_time)
%     frame_time=(0:5*60)*0.05;
%     for i=1:length(frame_time)
%     if rem(i,length(Im_new))
%         im_indx{i}=Im_new{rem(i,length(Im_new))}/255;
%     else
%         im_indx{i}=Im_new{length(Im_new)}/255;
%
%     end
% end
% end
%%Run over all channels
%        for AC(j)=find(MCS_electrodes==14) %%Run o ver all channels
for j=1:length(AC)
    %   for j=48
    for k=1:dim{AC(j)}


        if ~isempty(spike_ind_t{AC(j),k})
            for i=1:round(length(spike_ind_t{AC(j),k}))
                %                     Delay=frame_time(1);
                Delay=0;

                Event_times{AC(j),k}(i)=(spike_ind_t{AC(j),k}(i)/fs)-Delay;
                count_frame{AC(j)}(m)=1;

                for m=1:10
                    %              ind{m}=find(frame_time<=(Event_times{j}(i)-0.02-0.02*(m-1)) & frame_time>=(Event_times{j}(i)-0.04-0.02*(m-1)));
                    ind{m}=find(frame_time<=(Event_times{AC(j),k}(i)-(0.00+(FrameDuration*10^(-3)*(m-1))))); % Find Frames at least 40ms before spike.

                    %             ind_70=find(frame_time<=(Event_times{j}(i)-0.07) & frame_time>=(Event_times{j}(i)-0.120));
                    %             ind_120=find(frame_time<=(Event_times{j}(i)-0.1200));
                    %
                    if ~isempty( ind{m})
                        %                     for k=1:10
                        %                         if length(ind{m})-k>0
                        Frame{AC(j),k}(:,:,m)=(Frame{AC(j),k}(:,:,m)+double(im_indx{ind{m}(end)}));
                        count_frame{AC(j)}(m)= count_frame{AC(j)}(m)+1;
                        %                     end
                        %                 end
                    end

                end
            end
        end
        %
        max_val(AC(j),k)=max(max(max(Frame{AC(j),k}(:,:,:)))); % Calculate maximal value of each cluster for all 10 frames.
        std_val1(AC(j),k) = std(Frame{AC(j),k}(:,:,:)/255,1,'all');
    end
end

% range=max(max_val);
%display('Choose file path & name for analysed data')
%[file,path]=uiputfile('.mat','Choose file path & name for analysed data');

%save([path file(1:end-4)],'Frame','Aligned_Spikes','max_val','spike_ind_t','-v7.3')

%%
% for i = 1:size(Frame,1)
%     for k = 1:size(Frame,2)
% A{i,k} = Frame{i,k}/255; % Normlize to Pixel Amplitude
% A{i,k} = A{i,k}/size(spike_ind_t{i,k},2); % Normlize to spikes in each cluster
%     end
% end
% %% Spike-Triggered Covariance
% M = zeros(size(im_indx{1}(50:end,50:end),1),size(im_indx{1}(50:end,50:end),2));
% for i=1:length(im_indx)
% M = M+double(im_indx{i}(50:end,50:end));
% end
% M = M/i;
% Mcor = xcorr2(M); % Calculate the autocorrelation Mat of the prior (non spike related) stimiuli.
% Mcorr=Mcor(266:end-266,286:end-286);
% for j = 1:size(Frame,1)
%     for k = 1:size(Frame,2)
%         C{j,k} = zeros(size(Frame{1,1}(50:end,50:end),1),size(Frame{1,1}(50:end,50:end),1),10);
%         for i=1:length(spike_ind_t{j,k})
%               if spike_ind_t{j,k}(i)~=0
%                 Event_times{j,k}(i)=(spike_ind_t{j,k}(i)/fs);
%                 for m=1:size(Frame{j,k},3)
%                     ind{m}=find(frame_time<=(Event_times{j,k}(i)-0.04-0.04*(m-1)));
%                     if ~isempty( ind{m})
%                         S = double(im_indx{ind{m}(end)}(50:end,50:end));
%                         %C{j,k}(:,:,m) = C{j,k}(:,:,m) +  (S-Frame{j,k}(50:end,50:end,m))*((S-Frame{j,k}(50:end,50:end,m))');
%                         C{j,k}(:,:,m) = C{j,k}(:,:,m) + S*S';
%                     end
%                 end
%             end
%         end
%         C{j,k} = C{j,k}/(length(spike_ind_t{j,k}) - 1) - Mcorr;
%         for m=1:size(Frame{j,k},3)
%             [E{j,k}{m},D{j,k}{m}] = eig(C{j,k}(:,:,m));
%             [d,Ind] = sort(diag(D{j,k}{m}),'descend');
%             Ds{j,k}{m} = D{j,k}{m}(Ind,Ind);
%             Vs{j,k}{m} = E{j,k}{m}(:,Ind);
%         end
%     [Significance{j,k}] = NestedBootForSTC(im_indx,spike_ind_t{j,k},frame_time,Mcorr,Ds{j,k},fs);
%     end
% end
%% plot the maps 20msec
f=fspecial('gaussian',10,10);
% f=fspecial('gaussian',100,100 );
indx_elec=[];
% Channels_Psth= figure('name','Receptive Field');
% count=1;
% range = max(max(max(cell2mat(Frame))));
for i=1:length(AC)
    %if ~isempty(Aligned_Spikes{AC(i)})
        figure('Name',['Elec:',num2str(AC(i))]);
        count = 1;
        for k=1:dim{AC(i)}
            for m=1:10
                indx_elec =AC(i);

                subplot(3,10,count)
                temp=conv2(Frame{indx_elec(end),k}(50:end,90:end,m)/255,f,'same');
                temp1=Frame{indx_elec(end),k}(50:end,90:end,m);
                
                tempScale = temp/(max_val(indx_elec,k));
                % imagesc(temp/max(max(temp))-mean(mean(temp/max(max(temp))))
                %imagesc(temp,[0 (max_val(indx_elec,k)+0.001)/255]) % Added a minimal value to avoid [0 0] range in empty clusters
                imagesc(tempScale,[0 max(max(tempScale))])
                colormap(parula(1000))
                title([' Cluster:',num2str(k),' Frame:',num2str(m)])
                count=count+ 1;
            end
        end
    %end
end
%% Save Frame file for STA
frame = {Frame{12,1}(:,:,:),Frame{10,1}(:,:,:)}; % Insert the Channels and dimensions relevant for further analysis
save('C:\Users\Itay\Desktop\Yossi Mandel Lab\Thesis\Data Files\RF\29.11.21RF.mat','frame','AC');

%% calculate and plot STA
% sta= sta_Calc(Frame);
% [sta_cent, sta_surround]= sta_Calc(Frame);
% t_sta=0.4-(0:size(sta_surround,2)-1)/0.04;
% indx_elec=[];
%     count=16;
%
%     figure;
% while count>-1
%     for i=64:-1:1
%         subplot(8,8,64-i+1)
%         if rem(i,64)==1||rem(i,64)==8||rem(i,64)==57||rem(i,64)==0
%             cla(gca)
%             axis off
%         else
%             %         for k=1:length(chan)
%             if MEA_type==1
%                 indx_elec =[indx_elec;find(MCS_electrodes2(count)==MCS_electrodes)];
%
%                 subplot(8,8,i)
%                      plot(t_sta,(sta_surround(indx_elec(end),:)-mean(sta_surround(indx_elec(end),:)))/max((sta_surround(indx_elec(end),:)-mean(sta_surround(indx_elec(end),:)))),'b')
%                 hold on
%      plot(t_sta,(sta_cent(indx_elec(end),:)-mean(sta_cent(indx_elec(end),:)))/max((sta_cent(indx_elec(end),:)-mean(sta_cent(indx_elec(end),:)))),'r')
%
%                 ylabel('Amplitude[\muV]')
%                 xlabel('Time[msec]')
%                 %                                 title(['Elec',' ',num2str(MCS_electrodes2(count)),' ' ,'Z=',num2str(Imp(2,count))])
%                 %
%                 title(['Elec',' ',num2str(MCS_electrodes2(count))])
%                 count=count-1;
%                 %                 set(gca,'ButtonDownFcn', @mouseclick_callback5)
%                 %                 set(get(gca,'Children'),'ButtonDownFcn', @mouseclick_callback5)
%                 %
%                 set(gca,'ButtonDownFcn', @mouseclick_callback_sta)
%                 set(get(gca,'Children'),'ButtonDownFcn',@mouseclick_callback_sta)
%
%             else if MEA_type==2
%                     indx_elec =[indx_elec;find(MCS_electrodes2(count)==MCS_electrodes)];
%                     subplot(8,8,i)
%                     imagesc(Rec_50msec{indx_elec(end)})
%                     %                 hold on
%                     %                 plot(t,chan{indx_elec(end)},'r')
%                     %
%                     ylabel('Amplitude[\muV]')
%                     xlabel('Time[msec]')
%                     %                 title(['Elec','',num2str(Hybrid_electrodes(indx_elec(end))), ' Z=',num2str(Imp(indx_elec(end),2))], 'fontsize',8)
%                     title(['Elec','',num2str(Hybrid_electrodes(indx_elec(end)))])
%
%                     count=count-1;
%                     set(gca,'ButtonDownFcn', @mouseclick_callback_sta)
%                     set(get(gca,'Children'),'ButtonDownFcn', @mouseclick_callback_sta)
%                 end
%             end
%         end
%     end
% end
%
%
