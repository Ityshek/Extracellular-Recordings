function [SNR] = SNRCalc(Spike,raw_data,threshold,indx_spike)
Signal = [];
for i=1:length(Spike)
    Signal = [Signal;rms(Spike{i})];
end

window = floor(length(Spike{1})/2);
for c = 1:length(indx_spike)
raw_data(indx_spike(c)-window:indx_spike(c)+window) = NaN;
end

Noise = rms(rmmissing(raw_data));
SNR = 20*log10(mean(Signal)/Noise);
end