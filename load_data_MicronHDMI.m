function [raw_data, sampling_freq,stim_Data_seg,stim_sampling_rate,Begin_record,channelflag,stimulus_times, stimulus_indexes,answer,fname,pathname,SignalFiles,Dir]= load_data_MicronHDMI(c,answer,flag,pathname,SignalFiles,Dir)
if isempty(flag)
Dir=uigetdir('*.mat','Select a Folder to Load Raw Data From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = Dir;
end
count = 1;
    name=num2str(c);
    for i=1:length(SignalFiles)
        fname{i} = SignalFiles(i).name;
        [startIndex,endIndex] = regexp(fname{i},'_[0123456789._]*CP');
        if length(strsplit(fname{i}(startIndex+1:endIndex-3),'_')) >= 1 
            Valsstr = strsplit(fname{i}(startIndex+1:endIndex-3),'_');
            for v = 1:length(Valsstr)
            CPDs(v) =  str2num(Valsstr{v});
            end
        end
%             if isnan(CPDs)
%                 [startIndex,endIndex] = regexp(fname{i},'_\d+FrameDur');
%                 CPDs = str2double(fname{i}(startIndex+1:endIndex-8));
%             end
            data{i}=(load([Dir '\' SignalFiles(i).name]));
            if c<=9
                var_CSPK=['CSPK_00',name];
                var_CInPort=['CInPort_00',name];
                var_TimeBegin=['CSPK_00',name,'_TimeBegin'];
                var_TimeEnd=['CSPK_00',name,'_TimeEnd'];
            else
                var_CSPK=['CSPK_0',name];
                var_CInPort=['CInPort_0',name];
                var_TimeBegin=['CSPK_0',name,'_TimeBegin'];
                var_TimeEnd=['CSPK_0',name,'_TimeEnd'];
            end
            if isfield(data{i},var_CSPK)
                channelflag =0;
                raw_data_seg{i}=double(getfield(data{i},var_CSPK))*1.9;% convert to microvolts
                Begin_record_seg{i}=getfield(data{i},var_TimeBegin);
                End_record{i}=getfield(data{i},var_TimeEnd);
                sampling_freq=getfield(data{i},[var_CSPK,'_KHz'])*1000;
                stim_Data_seg{i}=data{i}.CInPort_001;
                stim_sampling_rate=data{i}.CInPort_001_KHz*1000;
            else
                raw_data = [];
                channelflag =1;
            end

        end
        raw_data = []; stim_Data = [];
        if ~isempty(raw_data_seg)
            for i=1:length(SignalFiles)
                raw_data =[raw_data,raw_data_seg{i}];
                raw_data_length(i) = length(raw_data_seg{i});
            end
            stimulus_times = stim_Data_seg{1}(1,1:2:end)/stim_sampling_rate -Begin_record_seg{1};
            %stimulus_indexes = round((stimulus_times)*sampling_freq);

            for i=2:length(SignalFiles)
                %b = stim_Data_seg{i}(1,1:2:end)/stim_sampling_rate +stimulus_times(end) - Begin_record_seg{i};
                b = (stim_Data_seg{i}(1,1:2:end)/stim_sampling_rate) +(sum(raw_data_length(1:i-1))/sampling_freq) - Begin_record_seg{i};       
                stimulus_times = [stimulus_times,b];
            end
            stimulus_indexes=round((stimulus_times)*sampling_freq);
            Begin_record = Begin_record_seg{1};
        end
count = count +1;    
end