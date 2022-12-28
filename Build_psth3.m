function [Psth,binsize_sec,smoothed_y,SpikeCount]=Build_psth3(Rast,fs)
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
    
    % Calculate number of spikes per 100ms between 10-110ms post stimulus
    temp2 = full(Rast(:,2*fs*0.01:12*fs*0.01));
    SpikeCount = sum(sum(temp2))/size(temp2,1);
    
    winsize = 4; % window size. 
    win = ones(winsize,1);
    smoothed_y = conv(Psth,win,'same')/winsize;