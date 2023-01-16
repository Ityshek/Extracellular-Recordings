clear all
 close all
clc
[filename, pathname]=uigetfile('*.mat');
load([pathname filename])
% for i=1:length(frame)
% frame{i} = frame{i}/max(max(max(frame{i})));
% end
Frame = frame';
[fwhmY, fwhmX, data_fit,data_fitx,Frame_int,CenterMassY, CenterMassX,F,data,datax]= Calc_Guassian_FitSingleunit(Frame,AC);

[e_center,e_surr1,e_surr2]=sta_center_SurroundSingleunit(fwhmY, fwhmX,Frame_int,CenterMassY, CenterMassX,Frame,data_fit,data_fitx,data,datax,AC);

[STA_cent,STA_Surr]=STA_Temporal_BehaviourSingleunit(e_center,e_surr1,e_surr2,Frame,AC);
