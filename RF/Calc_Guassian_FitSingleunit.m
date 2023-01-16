function [fwhmY, fwhmX, data_fit,data_fitx,Frame_int,CenterMassY, CenterMassX,F,data,datax]= Calc_Guassian_FitSingleunit(Frame,AC)

%find the frame with the highest std from the baseline for each channel

for j=1:length(Frame)
    for k=1:size(Frame,2)
        for m=1:10
            %  temp=conv2(Frame{j,k}(50:end,90:end,m),f,'same');
            temp=Frame{j,k}(50:end,90:end,m);
            
            %  temp2=temp/max(max(temp))-mean(mean(temp/max(max(temp))));
            std_base(j,m,k)=std(reshape(temp,size(temp,1)*size(temp,2),1),1);
            %         sta_cent(j,m)=max(max(std(Frame{j}(:,:,m))));
            %         sta_surround(j,m)=min(min(std(Frame{j}(:,:,m))));
        end
        [val indx]=max(std_base(j,:,k));
        %     Frame_int{j,k}=(Frame{j,k}(:,:,indx))-mean(mean(Frame{j,k}(:,:,indx)));
        Frame_int{j,k}=(Frame{j,k}(:,:,indx));
        
        val2(j,k)=max(max(Frame_int{j,k}));
        
    end
    
end
range=max(max(val2));
%% fit a guassian for each channel.
% for j=1:length(Frame)
%     for k=1:size(Frame,2)
%     frame_rec{j,k}=zeros(size(Frame_int{j,k}));
%     frame_vec{j,k}=reshape(Frame_int{j,k}(50:end,90:end),size(Frame_int{j,k}(50:end,90:end),1)*size(Frame_int{j,k}(50:end,90:end),2),1);
%     frame_vec2{j,k}=zeros(length(frame_vec{j,k}),1);
%     frame_rec_ind{j,k}= find(frame_vec{j,k}>=std( frame_vec{j,k}));
%      frame_vec2{j,k}(frame_rec_ind{j,k})=(frame_vec{j,k}(frame_rec_ind{j,k}));
%     frame_rec{j,k}=reshape(frame_vec2{j,k},size(Frame_int{j,k}(50:end,90:end),1),size(Frame_int{j,k}(50:end,90:end),2));
%     end
% end
%% preparing grid for the 2D guassian
% [X,Y] = meshgrid(1:size(Frame_int{j,k}(50:end,90:end),2),	1:size(Frame_int{j,k}(50:end,90:end),1));
[X,Y] = meshgrid(1:size(Frame_int{j,k}(50:end,90:end),2),	1:size(Frame_int{j,k}(50:end,90:end),1));

%% Fit to a guassian
% first find the center of mass of each frame of interest.
for j=1:length(Frame)
    %   for j=58
    
    for k=1:size(Frame,2)
        %    for k=2
        
        M=Frame_int{j,k}(50:end,90:end);
        thresh=0.88*max(max(M));
        [row col]=find(M>thresh&M<max(max(M)));
        %
        M_bin=zeros(size(M(:,:,1)));
        [rc,cc] = ndgrid(1:size(M,1),1:size(M,2));
        %
        M_bin(find(M>thresh))=1;
        
        
        Mt = sum(M_bin(:));
        c1 = sum(M_bin(:).* rc(:)) / Mt;
        c2= sum(M_bin(:) .* cc(:)) / Mt;
        CenterMassY(j,k)=c1;
        CenterMassX(j,k)=c2;
        [fwhmY(j,k), fwhmX(j,k), data_fit{j,k},data_fitx{j,k},data{j,k},datax{j,k}]=Gauss_FitSingleunit(M,CenterMassX(j,k),CenterMassY(j,k));
        % bulid depicted guassian
        %    if fwhmY(j,k)<250||fwhmX(j,k)<250
        va_exp=((X-CenterMassX(j,k)).^2)/(2*(fwhmX(j,k)/2.35).^2)+((Y-CenterMassY(j,k)).^2)/(2*(fwhmY(j,k)/2.35).^2);
        
        F{j,k}=exp(-va_exp);
        %    else
        %         F{j,k}=nan(size(Frame_int{j,k},1),size(Frame_int{j,k},2));
        %    end
    end
end

%% plot gaussian fits


f=fspecial('gaussian',10,10 );
% f=fspecial('gaussian',100,100 );

indx_elec=[];
% count=1;
% range = max(max(max(cell2mat(Frame))));
%ChannelPosition = [1,9,5,6,14,13,2,10,15,7,12,11,3,4,16,8];
%AC = [1:16];
ChannelPosition = [1,2];
% figure('Name',['Elec:',num2str(i)]);
for k=1:size(Frame,2)
    figure; 
    for i=1:length(ChannelPosition)
       
        
        
        
        count = ChannelPosition(i);
        indx_elec =i;
        
        subplot(1,2,count)
        
        
        imagesc( F{indx_elec(end),k})
        title(['Channel: ',num2str(AC(i))]);
        %
        %             temp=conv2(Frame{indx_elec(end),k}(50:end,90:end,m),f,'same');
        %             temp1=Frame{indx_elec(end),k}(50:end,90:end,m);
        %
        %            % imagesc(temp/max(max(temp))-mean(mean(temp/max(max(temp))))
        %
        %              imagesc(temp,[0 max_val(indx_elec,k)])
        
        %             ylabel('Amplitude[\muV]')
        %             xlabel('Time[msec]')
        %             title([' Cluster:',num2str(k),' Frame:',num2str(m)])
        
    end
end












