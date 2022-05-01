function Psth=Build_psth(Rast,fs)
 count1=1;
 
    binsize=10*10^-3*fs; 
    Psth=zeros(round(size(Rast,2)/binsize),1);
    binsize_sec=binsize/fs;
    for k=1:binsize:size(Rast,2)
        if k+binsize<size(Rast,2)
            temp=full(Rast(:,k:k+binsize-1));
            Psth(count1)=(((sum(sum(temp))))/size(Rast,1))/binsize_sec;
            count1=count1+1;
        end
    end

