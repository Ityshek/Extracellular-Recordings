function [CorrPlot] =  CorrFunc(t,ClusteredSpikeIdx,dim,sampling_freq)
if dim == 2
SpikeTrainMat = zeros(2,length(t));
for i=1:dim
for k = 1:length(t)
if find(ClusteredSpikeIdx{i} == k) ~= 0
    SpikeTrainMat(i,k) = 1;
end
end
DeltaAuto = 0.3; % Offset lag for autocorrelation 

AutoCorr(i,:) = xcorr(SpikeTrainMat(i,:),DeltaAuto*sampling_freq,'none');
t_AutoCorr = linspace(-DeltaAuto,DeltaAuto,sampling_freq*2*DeltaAuto+1);
SumSpikes = max(AutoCorr(i,:));
AutoCorr(i,find(AutoCorr(i,:) == max(AutoCorr(i,:)))) = 0;

end
DeltaCross = 0.3; % Offset lag for crosscorrelation

CrossCorr(1,:) = xcorr(SpikeTrainMat(1,:),SpikeTrainMat(2,:),DeltaCross*sampling_freq,'none');
SumSpikes1 = sum(CrossCorr(1,:));

CrossCorr(2,:) = xcorr(SpikeTrainMat(2,:),SpikeTrainMat(1,:),DeltaCross*sampling_freq,'none');
SumSpikes2 = sum(CrossCorr(2,:));
t_CrossCorr = linspace(-DeltaCross,DeltaCross,sampling_freq*2*DeltaCross+1);

CorrPlot = figure('Name','Correlation Plots');
for i = 1:dim
subplot(2,2,i)
bin_size = 0.005; % bin size in Sec
nbins = round((length(t_AutoCorr)/sampling_freq)/bin_size);
%histogram(AutoCorr(i,:),nbins)
plot(t_AutoCorr,AutoCorr(i,:))
xlabel('Offset [Sec]')
ylabel('Count')
title(['Autocorrelation Unit ',num2str(i)]);
subplot(2,2,i+2)
plot(t_CrossCorr,CrossCorr(i,:))
xlabel('Offset [Sec]')
ylabel('Count')
title(['Crosscorrelation refrence Unit ',num2str(i)]);
end

end

end