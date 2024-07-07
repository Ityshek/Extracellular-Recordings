function [raw,sampling_freq,stim_Data,stim_sampling_rate,Begin_record,Stim_Times,stimulus_indexes,StimIdx,StimFreqs,StimFreq,Check] =load_data_ConcateCFFNIR(ActiveChannels,Dir)
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = Dir;
count = 1;
[a,b] = regexp(SignalFiles(1).name,'_[0123456789_]*Hz');
valstr = strsplit(SignalFiles(1).name(a+1:b-2),'_');
for k=1:length(valstr)
    StimFreq(k) = str2num(valstr{k});
end
for c=ActiveChannels(count):ActiveChannels(end)
    name=num2str(c);
    for i=1:length(SignalFiles)
        fname{i} = SignalFiles(i).name;
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
            channelflag(c) =0;
            raw_data_seg{c}{i}=double(getfield(data{i},var_CSPK))*1.9;% convert to microvolts
            Begin_record_seg{c}{i}=getfield(data{i},var_TimeBegin);
            End_record{c}{i}=getfield(data{i},var_TimeEnd);
            sampling_freq=data{i}.CSPK_001_KHz*1000;
            stim_Data_seg{c}{i}=data{i}.CAI_001;
            stim_sampling_rate=data{i}.CAI_001_KHz*1000;
        else
            raw_data{c} = [];
            channelflag(c) =1;
        end

    end
    raw_data{c} = []; stim_Data{c} = [];
    if ~isempty(raw_data_seg{c})
        for i=1:length(SignalFiles)
            raw_data{c} =[raw_data{c},raw_data_seg{c}{i}];
            raw_data_length(i) = length(raw_data_seg{c}{i});
        end

        Begin_record{c} = Begin_record_seg{c}{1};
    end
    Cat_Stim_Seg{c} = [];
    % Concat Trigger Channel
    for i=1:length(SignalFiles)
        Cat_Stim_Seg{c} = [Cat_Stim_Seg{c},stim_Data_seg{c}{i}];
    end
    % Align Trigger Channel for detection
    for i=1:length(Cat_Stim_Seg{c})
        if Cat_Stim_Seg{c}(i) <10 && StimChannel(i) > -10
            Cat_Stim_Seg{c}(i) = NaN;
        end
        if Cat_Stim_Seg{c}(i) ~= 0
            Cat_Stim_Seg{c}(i) = Cat_Stim_Seg{c}(i)-max(Cat_Stim_Seg{c});
        end
    end
    %% Retrive Triggers from TTL
    % for i = 1:length(data)
    %     sampling_rate_stim_led=1000*data{i}.CTTL_002_KHz;
    %     %converting trigger indecies to times in secs
    %     stimulus_times_Up{i}=(data{i}.CTTL_002_Up(1,:))/(sampling_rate_stim_led)+data{i}.CTTL_002_TimeBegin-data{i}.CLFP_001_TimeBegin;
    %     stimuus_indexes_Up{i}=stimulus_times_Up{i}*sampling_freq;
    %
    %     stim_index{i}=stimuus_indexes_Up{i};
    %     Stimulus_times{i}=stim_index{i}/sampling_freq;
    % end
    % stimulus_times = [Stimulus_times{1},Stimulus_times{2}+Stimulus_times{1}(end)+(data{1, 2}.CLFP_001_TimeBegin-data{1, 1}.CTTL_002_TimeEnd)];
    % stimulus_indexes = stimulus_times*sampling_freq;

    %%
    stim_thresh = -1800;
    stimulus_times=(find(circshift(Cat_Stim_Seg{c},[0 1])<stim_thresh&Cat_Stim_Seg{c}>stim_thresh))/stim_sampling_rate;
    FirstPulse = StimFreq(1);
    [a,b] = regexp(SignalFiles(1).name,'_[0123456789]*ms');
    PulseWidth = str2num(SignalFiles(1).name(a+1:b-2));
    IPI = (1000/FirstPulse - PulseWidth)/1000;
    stimulus_times = [(stimulus_times(1)-IPI),stimulus_times];
    stimulus_times = stimulus_times(1:end-1);
    %
    % Check if stimulus times marks the correct peaks in the trigger
    figure(); plot(Cat_Stim_Seg{c}); hold on;
    even=nan(length(Cat_Stim_Seg{c}),1);
    even(round(stimulus_times*stim_sampling_rate))=-2000; plot(even,'*r')

    Check = stimulus_times;

    NumReps = 40;
    CountStim=1;CountStim2=StimFreq(1)+1;
    StimFreqs = StimFreq;
    for i=1:NumReps-1
        StimFreqs = [StimFreqs,StimFreq];
    end
    if length(stimulus_times) == 5080
        StimIdx = []; StimIdx(1) = 1; 
        Stim_Times(1) = stimulus_times(1);
        for k=2:length(StimFreqs)
            StimIdx(k) = StimIdx(end)+StimFreqs(k-1);
            Stim_Times(end+1) = stimulus_times(StimIdx(k));
        end
    else
        StimFreqs = StimFreqs(2:end);
        StimIdx = []; StimIdx(1) = 2; StimIdx(2) = StimFreq(1)+StimIdx(1)+1;
        Stim_Times(1) = stimulus_times(StimFreq(1));
        for k=3:length(StimFreqs)
            StimIdx(k) = StimIdx(end)+StimFreqs(k-1);
            Stim_Times(end+1) = stimulus_times(StimIdx(k));
        end
    end
    hold on; plot(Stim_Times*2750,ones(1,length(Stim_Times))*(-1800),'*g')
    stimulus_indexes = round(Stim_Times*sampling_freq);
    count = count +1;
    %Stim_Times = stimulus_times; stimulus_indexes = int32(stimulus_indexes); StimIdx = int32(stimulus_indexes);
    %hold on; plot(Stim_Times*stim_sampling_rate,ones(1,278)*1500,'g*') % Check final stim locations
end

raw = raw_data{ActiveChannels};
end