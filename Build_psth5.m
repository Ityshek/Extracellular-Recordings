function [Psth,binsize_sec]=Build_psth5(Rast,fs,BinSize)
 count1=1;
 
    binsize=BinSize*fs;
    WindowSize = size(Rast,2);
    Psth=zeros(round(WindowSize/binsize),1);
    binsize_sec=binsize/fs;
    for k=1:binsize:WindowSize
        if k+binsize-1<=WindowSize
            temp=full(Rast(:,k:k+binsize-1));
            Psth(count1)=(((sum(sum(temp))))/size(Rast,1))/binsize_sec;
            count1=count1+1;
        end
    end

