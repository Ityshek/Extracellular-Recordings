function [raw_data,sampling_freq,stim_Data,stim_sampling_rate,Begin_record,stimulus_times,stimulus_indexes,CPDs] =load_data_ConcateMultiUnit(ActiveChannels)
Dir=uigetdir('*.mat','Select a Folder to Load Raw Data From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = Dir;
for c=1:length(ActiveChannels)% 16 Channels
    name=num2str(ActiveChannels(c));
    for i=1:length(SignalFiles)
        fname{i} = SignalFiles(i).name;
        [startIndex,endIndex] = regexp(fname{i},'_\d+CPD');
        if length(str2mat((fname{i}(startIndex+1:endIndex-3))))<4
            CPDs(i) = str2double(fname{i}(startIndex+1:endIndex-3))*0.01;
        else
            CPDs(i) = str2double(fname{i}(startIndex+1:endIndex-3))*0.001;
        end
            if isempty(CPDs)
                [startIndex,endIndex] = regexp(fname{i},'_\d+FrameDur');
                CPDs = str2double(fname{i}(startIndex+1:endIndex-3));
            end
            data{i}=(load([Dir '\' SignalFiles(i).name]));
            if ActiveChannels(c)<=9
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
                channelflag(ActiveChannels(c)) =0;
                raw_data_seg{ActiveChannels(c)}{i}=double(getfield(data{i},var_CSPK))*1.9;% convert to microvolts
                Begin_record_seg{ActiveChannels(c)}{i}=getfield(data{i},var_TimeBegin);
                End_record{ActiveChannels(c)}{i}=getfield(data{i},var_TimeEnd);
                sampling_freq=data{i}.CSPK_001_KHz*1000;
                stim_Data_seg{ActiveChannels(c)}{i}=data{i}.CInPort_001;
                stim_sampling_rate=data{i}.CInPort_001_KHz*1000;
            else
                raw_data{ActiveChannels(c)} = [];
                channelflag(ActiveChannels(c)) =1;
            end

        end
        raw_data{ActiveChannels(c)} = []; stim_Data{ActiveChannels(c)} = [];
        if ~isempty(raw_data_seg{ActiveChannels(c)})
            for i=1:length(SignalFiles)
                raw_data{ActiveChannels(c)} =[raw_data{ActiveChannels(c)},raw_data_seg{ActiveChannels(c)}{i}];
                raw_data_length(i) = length(raw_data_seg{ActiveChannels(c)}{i});
            end
            stimulus_times{ActiveChannels(c)} = stim_Data_seg{ActiveChannels(c)}{1}(1,1:2:end)/stim_sampling_rate -Begin_record_seg{ActiveChannels(c)}{1};
            stimulus_indexes{ActiveChannels(c)} = round((stimulus_times{ActiveChannels(c)})*sampling_freq);

            for i=2:length(SignalFiles)
                b = stim_Data_seg{ActiveChannels(c)}{i}(1,1:2:end)/stim_sampling_rate +stimulus_times{ActiveChannels(c)}(end) - Begin_record_seg{ActiveChannels(c)}{i};
                stimulus_times{ActiveChannels(c)} = [stimulus_times{ActiveChannels(c)},b];
            end
            stimulus_indexes{ActiveChannels(c)}=round((stimulus_times{ActiveChannels(c)})*sampling_freq);
            Begin_record{ActiveChannels(c)} = Begin_record_seg{ActiveChannels(c)}{1};
        end
    end
end