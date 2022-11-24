function [Psth,binsize_sec,smoothed_y]=Build_psth3(Rast,fs)
 count1=1;
 
    binsize=5*10^-3*fs; 
    Psth=zeros(round(size(Rast,2)/binsize),1);
    binsize_sec=binsize/fs;
    for k=1:binsize:size(Rast,2)
        if k+binsize-1<=size(Rast,2)
            temp=full(Rast(:,k:k+binsize-1));
            Psth(count1)=(((sum(sum(temp))))/size(Rast,1))/binsize_sec;
            count1=count1+1;
        end
    end
    winsize = 4; % window size. 
    win = ones(winsize,1);
    smoothed_y = conv(Psth,win,'same')/winsize;