function [raw_data, sampling_freq,stim_Data,stim_sampling_rate,Begin_record,channelflag] =load_data_MultiUnit(c,fname,pathname)

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
sampling_freq=data.CSPK_001_KHz*1000;
stim_Data=data.CInPort_001;
stim_sampling_rate=data.CInPort_001_KHz*1000;
Begin_record=data.CSPK_001_TimeBegin;