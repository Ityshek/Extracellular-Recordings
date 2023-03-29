function [Psth,binsize_sec,smoothed_y,SpikeCount,Label]=Build_psth3(Rast,fs,CountWindow)
prompt = {'Choose Spike Calc: (1 = Hz | 2 = Count)','Choose Window Size [ms]:'};
definpt = {'2','10'}; dlgtitle = 'Input'; dims = [1 35];
Flag = inputdlg(prompt,dlgtitle,dims,definpt);
if Flag{1} == '1'
% Calc PSTH in Hz
count1=1;
binsize=str2num(Flag{2})*10^-3*fs;
%Psth=zeros(round(size(Rast,2)/binsize),1);
binsize_sec=binsize/fs;
for k=1:binsize:size(Rast,2)
    if k+binsize-1<=size(Rast,2)
        temp=full(Rast(:,k:k+binsize-1));
        Psth(count1)=(((sum(sum(temp))))/size(Rast,1))/binsize_sec;
        count1=count1+1;
    end
end
Psth = circshift(Psth,1);
winsize = 5; % window size.
win = ones(winsize,1);
smoothed_y = conv(Psth,win,'same')/winsize;
Label = 'Firing Rate [Hz]';
elseif Flag{1} == '2'
% Calac PSTH in num of spikes in bin
count1=1;
binsize=str2num(Flag{2})*10^-3*fs;
%Psth2=zeros(round(size(Rast,2)/binsize),1);
binsize_sec=binsize/fs;
for k=1:binsize:size(Rast,2)
    if k+binsize-1<=size(Rast,2)
        temp=full(Rast(:,k:k+binsize-1));
        Psth(count1)=(((sum(sum(temp))))/size(Rast,1));
        count1=count1+1;
    end
end

Psth = circshift(Psth,1);
winsize = 5; % window size.
win = ones(winsize,1);
smoothed_y = conv(Psth,win,'same')/winsize;
Label = ['Spike Count [',num2str(Flag{2}),'ms]'];
end
% Calculate number of spikes per 200ms between 10-210ms post stimulus
temp2 = full(Rast(:,CountWindow(1)*fs*0.01:CountWindow(2)*fs*0.01));
SpikeCount = sum(sum(temp2))/size(temp2,1);