function [Psth,binsize_sec,smoothed_y]=Build_psth4(Rast,fs,CountWindow,Flag)
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

%Psth = circshift(Psth,1);
winsize = 2; % window size.
win = ones(winsize,1);
smoothed_y = conv(Psth,win,'same')/winsize;
%smoothed_y = conv(Psth,win,'same');

