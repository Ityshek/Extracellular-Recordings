function [SNR] = SNRCalc(Spike,raw_data)
Signal = [];
for i=1:length(Spike)
    Signal = [Signal;rms(Spike{i})];
        end
Noise = rms(raw_data);
SNR = 20*log10(mean(Signal)/Noise);
end