function [e_center,e_surr1,e_surr2]=sta_center_SurroundSingleunit(fwhmY, fwhmX,Frame_int,CenterMassY, CenterMassX,Frame,data_fit,data_fitx,data,datax,AC)

%% for center , find mean of pixels in a radius of 1 sigma
[X,Y] = meshgrid(1:size(Frame_int{1,1}(50:end,90:end),2),	1:size(Frame_int{1,1}(50:end,90:end),1));
global data_fit data_fitx data datax




f=fspecial('gaussian',10,10 );
% f=fspecial('gaussian',100,100 );
%ChannelPosition = [1,9,5,6,14,13,2,10,15,7,12,11,3,4,16,8];
ChannelPosition = [1,2];
%AC = [1:16];

%Channels_Psth= figure('name','Receptive Field');
% count=1;
% range = max(max(max(cell2mat(Frame))));

%     figure('Name',['Elec:',num2str(i)]);
  
     for k=1:size(Frame,2)


for i=1:length(ChannelPosition)
    figure
      count = ChannelPosition(i);
      indx_elec =i;

%             subplot(1,2,count)
            
            
                  imagesc( conv2(Frame_int{indx_elec(end),k}(50:end,90:end),f,'same'))
                  title(['Channel: ',num2str(AC(i))]);
                    hold on
                    e_center{indx_elec(end),k}=((X-CenterMassX(indx_elec(end),k))/(1*(fwhmX(indx_elec(end),k)/2.35))).^2+((Y-CenterMassY(indx_elec(end),k))/(1*(fwhmY(indx_elec(end),k)/2.35))).^2<=1;
                    e_surr1{indx_elec(end),k}=((X-CenterMassX(indx_elec(end),k))/(1.5*(fwhmX(indx_elec(end),k)/2.35))).^2+((Y-CenterMassY(indx_elec(end),k))/(1.5*(fwhmY(indx_elec(end),k)/2.35))).^2<=1;
                    e_surr2{indx_elec(end),k}=((X-CenterMassX(indx_elec(end),k))/(2*(fwhmX(indx_elec(end),k)/2.35))).^2+((Y-CenterMassY(indx_elec(end),k))/(2*(fwhmY(indx_elec(end),k)/2.35))).^2<=1;
                    
                    hold on
                    %contour(e_center{indx_elec(end),k}+e_surr1{indx_elec(end),k}+e_surr2{indx_elec(end),k})
                    contour(e_center{indx_elec(end),k}+e_surr1{indx_elec(end),k})                   
      set(gca,'ButtonDownFcn',@mouseclick_cross_section_singleunit)
                    set(get(gca,'Children'),'ButtonDownFcn',@mouseclick_cross_section_section_singleunit)
     end
end








