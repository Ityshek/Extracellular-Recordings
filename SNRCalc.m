function [SNR] = SNRCalc(Spike,raw_data,threshold)
Signal = [];
for i=1:length(Spike)
    Signal = [Signal;rms(Spike{i})];
        end
Noise = rms(raw_data(find(raw_data>threshold)));
SNR = 20*log10(mean(Signal)/Noise);
end