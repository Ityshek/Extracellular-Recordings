 function mouseclick_cross_section_singleunit(gcbo,eventdata)
      % the arguments are not important here, they are simply required for
      % a callback function. we don't even use them in the function,
      % but Matlab will provide them to our function, we we have to
      % include them.
      ChannelPosition = [1,9,5,6,14,13,2,10,15,7,12,11,3,4,16,8];

      global data_fit data_fitx data datax % first we get the point that was clicked on
      cP = get(gca,'Currentpoint');
      x = cP(1,1);
      y = cP(1,2);
       H=get(gca,'title');
% 
      range=caxis;
%   figure; 
      % Now we find out which mouse button was clicked, and whether a
      % keyboard modifier was used, e.g. shift or ctrl
      switch get(gcf,'SelectionType')
          case 'open'   % Double-click any mouse button.
            
              h=get(gcbo,'Children');
                cb = getappdata(gca,'ColorbarPeerHandle');
%               x=get(h,'xdata');
figure
%               y=get(h,'ydata');



subplot(3,3,[4,5,7,8])

elec_indx=find(str2num(H.String(end-1:end))==ChannelPosition);
imagesc(h(2).CData,range)
hold on
contour(h(1).ZData)
 title(H.String)
subplot(3,3,[1:2])
plot([1:length(data_fitx{elec_indx,2})],data_fitx{elec_indx,2},'-b',[1:length(datax{elec_indx,2})],datax{elec_indx,2}/max(datax{elec_indx,2})-mean(datax{elec_indx,2}/max(datax{elec_indx,2})),'*r')

axis off
subplot(3,3,[6,9])

% plot(data_fit{elec_indx,2},[length(data_fit{elec_indx,2}):-1:1],'linewidth',2)
plot(data_fit{elec_indx,2},[length(data_fit{elec_indx,2}):-1:1],'-b',data{elec_indx,2}/max(data{elec_indx,2})-mean(data{elec_indx,2}/max(data{elec_indx,2})),[length(data_fit{elec_indx,2}):-1:1],'*r')

axis off



%               ylim([0 5000])
%                ylabel('Rate[Spike/sec]')
%             xlabel('Time[msec]')
                       

          end
      % get and set title handle
     
  end
