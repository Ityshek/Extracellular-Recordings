function [raw_data, sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag,stimulus_times, stimulus_indexes] = SUload_data_Micron(c,fname,pathname)
StimType = str2num(cell2mat(inputdlg('Insert 2 for Vis / 1 for NIR')));
name=num2str(c);
data=load([pathname fname]);
if c<=9
    var_str=['CSPK_00',name];
else
    var_str=['CSPK_0',name];
end
if isfield(data,var_str)
    channelflag =0;
    raw_data=double(getfield(data,var_str))*1.9;% convert to microvolts
else
    raw_data = [];
    channelflag =1;
end
if StimType == 2
    StimChannel = data.CAI_002;
    stim_thresh=1000;
else
    StimChannel = data.CAI_001;
    stim_thresh=500;
end
sampling_freq=data.CSPK_001_KHz*1000;
stim_Data=StimChannel; 
stim_sampling_rate=data.CAI_002_KHz*1000;
t_stim=(0:length(StimChannel)-1)/stim_sampling_rate;
stimulus_times=(find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate;
stimulus_indexes=round(stimulus_times*sampling_freq);
stim=zeros(length(t_stim),1);
stim(round(stimulus_times*stim_sampling_rate))=  stim_thresh;

Begin_record=data.CSPK_001_TimeBegin;
end