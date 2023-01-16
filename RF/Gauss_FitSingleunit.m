function [fwhmY, fwhmX, data_fit,data_fitx,data,datax]=Gauss_FitSingleunit(M,CenterX,CenterY)
% if round(CenterY-200)>0&round(CenterX-200)>0&round(CenterX+200)<size(M,2)&round(CenterY+200)<size(M,1)
%     m=M(round(CenterY-200):round(CenterY+200),round(CenterX-200):round(CenterX+200));
% else if round(CenterY-200)<0||round(CenterX-200)<0
%     m=M(round(CenterY):round(CenterY+200),round(CenterX):round(CenterX+200));
% else if round(CenterY+200)>size(M,1)||round(CenterX+300)>size(M,2)
%     m=M(round(CenterY-200):round(CenterY),round(CenterX-200):round(CenterX));
%     end
%     end
%     
% end
%  figure;
%    imagesc(m)
%    colormap('gray')
%   [x y]=ginput(1);
  data=M(1:end,round(CenterX));
%     plot(data,'.')
  peak_data=max(data);
  baseline=min(data);
  data_new=double(data-baseline);
  x_axis=1:length(data_new);
  Range=0.5*(peak_data-baseline);
 indicies=find(data_new>=Range-0.5&data_new<=Range+0.5);
  func2=@(x_var)sum((data_new'/max(data_new)-x_var(3)*((1/(2*pi*x_var(1)^2)^0.5).*exp(-(x_axis-x_var(2)).^2/(2*x_var(1)^2)))).^2);
  [x_var,eval]=fminsearch(func2,[  baseline,CenterX,  baseline]);
  data_fit=(x_var(3)*(1/(2*pi*x_var(1)^2)^0.5).*exp(-(x_axis-x_var(2)).^2/(2*x_var(1)^2)));
%   figure;
%  plot(x_axis,data_new/max(data_new),'*',x_axis,data_fit)
fwhmY=abs(x_var(1))*(8*log(2))^0.5;
%% Repeat for X-direction (orhtogonol )

 datax=M(round(CenterY),1:end);
%     plot(datax,'.')
  peak_datax=max(datax);
  baselinex=min(datax);
  data_newx=double(datax-baselinex);
  y_axis=1:length(data_newx);
  Rangex=0.5*(peak_datax-baseline);
 indiciesx=find(data_newx>=Rangex-0.5&data_newx<=Rangex+0.5);
  func3=@(x_var2)sum((data_newx/max(data_newx)-x_var2(3)*((1/(2*pi*x_var2(1)^2)^0.5).*exp(-(y_axis-x_var2(2)).^2/(2*x_var2(1)^2)))).^2);
  [x_var2,eval]=fminsearch(func3,[baselinex,CenterY,baselinex]);
  data_fitx=(x_var2(3)*(1/(2*pi*x_var2(1)^2)^0.5).*exp(-(y_axis-x_var2(2)).^2/(2*x_var2(1)^2)));
%   figure;
%   plot(y_axis,data_newx/max(data_newx),'*',y_axis,data_fitx)
fwhmX=abs(x_var2(1))*(8*log(2))^0.5;