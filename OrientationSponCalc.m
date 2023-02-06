function [spon] = OrientationSponCalc(spikeidx,stimulus_times,sampling_freq,BlankScreenSec)
spike_times = spikeidx/sampling_freq;
for i = 1:length(stimulus_times)
SpikeCount(i) = length(find(spike_times>stimulus_times(i)&spike_times<stimulus_times(i)+BlankScreenSec));

end

SponWindow = mean(diff(spike_times)) - BlankScreenSec;
spon = mean(SpikeCount);
end