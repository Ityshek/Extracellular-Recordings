function [raw_data, sampling_freq,StimChannel,stim_sampling_rate,Begin_record,channelflag,stimulus_times, stimulus_indexes,answer]= load_data_Micron(c,fname,pathname,answer,NumReps)

StimType = answer{1}; DistortTrigger = answer{2};  Type = answer{3};
name=num2str(c);
data=load([pathname fname]);
if isfield(data,'data')
    data = data.data;
end
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
%% Check for NIR condition
if StimType == '1'
    StimChannel = data.CAI_001;
    stim_sampling_rate=data.CAI_001_KHz*1000;
    sampling_freq=getfield(data,[var_str,'_KHz'])*1000;
    stim_Data=StimChannel;

    if Type == '1' % Check for Alternating Condition
        StimChannel = data.CAI_002;
        stim_thresh = 1000;
        NumFrames = 2;
    elseif Type == '2' % Check for CFF condition
        %         for i=1:length(StimChannel)
        %             if StimChannel(i) <10 && StimChannel(i) > -10
        %                 StimChannel(i) = NaN;
        %             end
        %             if StimChannel(i) ~= 0
        %                 StimChannel(i) = StimChannel(i)-max(StimChannel);
        %             end
        %         end
        [a,b] = regexp(fname,'_[0123456789]*Hz');
        NumFrames = str2num(fname(a+1:b-2));
        [a,b] = regexp(fname,'_[0123456789]*Reps');
        Reps = str2num(fname(a+1:b-4));
        %         for s=3:abs(min(StimChannel))
        %             NumTrigs = NumFrames*Reps; % Multiply by repetitions
        %             stim_thresh = -s;
        %             if length((find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate) == NumTrigs
        %                 break
        %             end
        %         end
    else
        stim_thresh = min(StimChannel)+50;
        NumFrames = 1;
        %stim_thresh = 5;
    end
    [startIndex,endIndex] = regexp(fname,'_\d*Hz');
    StimFreq = str2num(fname(startIndex+1:endIndex-2));
    [startIndex,endIndex] = regexp(fname,'_\d*ms');
    %     if StimFreq == 64
%     for s = min(StimChannel)+10: max(StimChannel)-10
%         stim_thresh = s;
%         stimulus_times=(find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate;
%         if length(stimulus_times) == StimFreq*NumReps
%             break
%         end
%     end
    %     else
    %         stim_thresh = min(StimChannel)+100;
    %         stimulus_times=(find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate;
    %     end
    %stimulus_times = stimulus_times(1:round(length(stimulus_times)/NumFrames,0))
    %if Type == '1'
    stimulus_times=(find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate;
    stimulus_times = stimulus_times(1:NumFrames:end);
    %else
    %stimulus_times = stimulus_times(1:NumFrames:end);
    %end
    stimulus_indexes=round(stimulus_times*sampling_freq);
    %stim=zeros(length(t_stim),1);
    %stim(round(stimulus_times*stim_sampling_rate))=  stim_thresh;
    %% Check for VIS condition
elseif StimType == '2'
        StimChannel = data.CAI_002;
    %     StimChannel = StimChannel - min(StimChannel);
    stim_thresh=30;
    if Type == '-' % Check for Null / Intensity Condition
        stim_thresh=max(StimChannel)-10; NumFrames = 2;
        sampling_freq=getfield(data,[var_str,'_KHz'])*1000;
        stim_sampling_rate=data.CAI_002_KHz*1000;
        t_stim=(0:length(StimChannel)-1)/stim_sampling_rate;
    
    for s = double(min(StimChannel)+10):0.1:double(max(StimChannel)-1)
        stim_thresh = s;
        stimulus_times=(find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate;
        stimulus_times = stimulus_times(1:NumFrames:end); A = round(diff(stimulus_times),2);
        if sum(~ismember(A,1)) == 0 && ~isempty(A) % Check for all 
            break
        end
    end
        stimulus_times=(find(circshift(StimChannel,[0 -1])<stim_thresh&StimChannel>stim_thresh))/stim_sampling_rate;
        stimulus_times = stimulus_times(1:NumFrames:end-1);
        % stimulus_times = stimulus_times(1+NumFrames:NumFrames:41*NumFrames); 
        stimulus_indexes=round(stimulus_times*sampling_freq);
    elseif Type == '2' % Check for CFF condition
%         [a,b] = regexp(fname,'_[0123456789]*Hz');
%         NumFrames = str2num(fname(a+1:b-2))*2;
%         sampling_freq=data.CSPK_017_KHz*1000;
%         stim_sampling_rate=data.CAI_002_KHz*1000;
% 
%         if NumFrames == 128
%             val = max(StimChannel); StimChannel = abs(StimChannel-val);
%             val2 = max(StimChannel);
%             for i=1:length(StimChannel)
%                 if  StimChannel(i) > val2 - 5
%                     stimchannel(i) = 0;
%                 else
%                     stimchannel(i) = 40;
%                 end
%             end
%             stim_thresh=20;
%             stimulus_time=(find(circshift(stimchannel,[0 1])<stim_thresh&stimchannel>stim_thresh))/stim_sampling_rate;
%             stimulus_times = linspace(stimulus_time,stimulus_time+((NumReps-1)*2),NumReps);
%             stimulus_indexes=round(stimulus_times*sampling_freq);
% 
%         else
%             %         [startIndex,endIndex] = regexp(fname,'_\d*Hz');
%             %         StimFreq = str2num(fname(startIndex+1:endIndex-2));
%             t_stim=(0:length(StimChannel)-1)/stim_sampling_rate;
%             for s = min(StimChannel)+10: max(StimChannel)-10
%                 stim_thresh = s;
%                 stimulus_times=(find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate;
%                 A = max(diff(stimulus_times(1:NumFrames:end))); B = min(diff(stimulus_times(1:NumFrames:end)));
%                 if A < 2.01 && B > 1.99
%                     break
%                 end
%             end
%             %stim_thresh = 1000;
%             stimulus_times=(find(circshift(StimChannel,[0 -1])<stim_thresh&StimChannel>stim_thresh))/stim_sampling_rate;
%             stimulus_times= stimulus_times(1:NumFrames:end);
%             %                 if NumFrames == 128
%             %                     stimulus_times = stimulus_times(1:NumFrames-3:end);
%             %                 else
%             %                     EndPoint = round(length(stimulus_times)/NumFrames)*NumFrames;
%             %                     stimulus_times = stimulus_times(1:NumFrames:EndPoint);
%             %                 end
%             stimulus_indexes=round(stimulus_times*sampling_freq);
%             stim=zeros(length(t_stim),1);
%             stim(round(stimulus_times*stim_sampling_rate))=  stim_thresh;
%         end
        
        sampling_freq=data.CSPK_017_KHz*1000;
        stim_sampling_rate=data.CAI_002_KHz*1000;
        [a,b] = regexp(fname,'_[0123456789]*Hz');
        NumFrames = str2num(fname(a+1:b-2));
        [a,b] = regexp(fname,'_[0123456789]*Reps');
        Reps = str2num(fname(a+1:b-4));
            [startIndex,endIndex] = regexp(fname,'_\d*Hz');
    StimFreq = str2num(fname(startIndex+1:endIndex-2));
    [startIndex,endIndex] = regexp(fname,'_\d*ms');
    %     if StimFreq == 64
    for s = min(StimChannel)+10: max(StimChannel)-10
        stim_thresh = s;
        stimulus_times=(find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate;
        stimulus_times = stimulus_times(1:NumFrames*2:end); A = round(diff(stimulus_times),1);
        if sum(~ismember(A,2)) == 0 % Check for all 
            break
        end
    end
    stimulus_times=(find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate;
    stimulus_times = stimulus_times(1:NumFrames*2:end);
    stimulus_indexes=round(stimulus_times*sampling_freq);
    end
    %% Check for Micron HDMI condition
elseif StimType == '3'
    sampling_freq=getfield(data,[var_str,'_KHz'])*1000; stim_sampling_rate=data.CInPort_001_KHz*1000;
    StimChannel = data.CInPort_001; Begin_record=getfield(data,[var_str,'_TimeBegin']);
    stimulus_times=((StimChannel(1,1:2:end))/stim_sampling_rate)-Begin_record;
    stimulus_indexes=round((stimulus_times)*sampling_freq);
    % else
    %     stim_thresh = 1000;
    %     NumFrames = 2;
    %     sampling_freq=getfield(data,[var_str,'_KHz'])*1000;
    %     stim_sampling_rate=data.CAI_002_KHz*1000;
    %     t_stim=(0:length(StimChannel)-1)/stim_sampling_rate;
    %     stimulus_times=(find(circshift(StimChannel,[0 -1])<stim_thresh&StimChannel>stim_thresh))/stim_sampling_rate;
    %     stimulus_times = stimulus_times(1:NumFrames:end);
    %     stimulus_indexes=round(stimulus_times*sampling_freq);
    %     stim=zeros(length(t_stim),1);
    %     stim(round(stimulus_times*stim_sampling_rate))=  stim_thresh;
    if DistortTrigger == '1'
        First_stim = min((find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate);
        Last_stim = max((find(circshift(StimChannel,[0 -1])>stim_thresh&StimChannel<stim_thresh))/stim_sampling_rate);
        [startIndex,endIndex] = regexp(fname,'_\d*Hz');
        Trigger_Freq = str2double(fname(startIndex+1:endIndex-2));
        stimulus_times = [First_stim:1/Trigger_Freq:Last_stim];
        stimulus_indexes=round(stimulus_times*sampling_freq);
        stim=zeros(length(t_stim),1);
        stim(round(stimulus_times*stim_sampling_rate))=  stim_thresh;
    end

end
Begin_record=getfield(data,[var_str,'_TimeBegin']);
%%
