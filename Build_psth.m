function [Psth,binsize_sec,smoothed_y,SpikeCount,Label,window]=Build_psth(Rast,fs)
%  count1=1;
%  
%     binsize=5*10^-3*fs; 
%     Psth=zeros(round(size(Rast,2)/binsize),1);
%     binsize_sec=binsize/fs;
%     for k=1:binsize:size(Rast,2)
%         if k+binsize<size(Rast,2)
%             temp=full(Rast(:,k:k+binsize-1));
%             Psth(count1)=(((sum(sum(temp))))/size(Rast,1))/binsize_sec;
%             count1=count1+1;
%         end
%     end
%     
%     % Calculate number of spikes per 150ms between 10-160ms post stimulus
%     temp2 = full(Rast(:,2*fs*0.01:17*fs*0.01));
%     SpikeCount = sum(sum(temp2))/size(temp2,1);
%     
%     
%     winsize = 4; % window size. 
%     win = ones(winsize,1);
%     smoothed_y = conv(Psth,win,'same')/winsize;
%     PSTHbinsize =binsize/fs*10^3;
    %%
prompt = {'Choose Spike Calc: (1 = Hz | 2 = Count)','Choose Bin Size [ms]:','choose Window in ms post-stimulus:' };
definpt = {'2','10','[10 210]'}; dlgtitle = 'Input'; dims = [1 35];
Flag = inputdlg(prompt,dlgtitle,dims,definpt);
window = str2num(Flag{3});
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
temp2 = full(Rast(:,window(1)*0.001*fs+440:window(2)*0.001*fs+440));
SpikeCount = sum(sum(temp2))/size(temp2,1);