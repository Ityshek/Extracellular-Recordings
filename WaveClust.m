function  [idx,col] = WaveClust(waveforms,Idxwaves)
% Define parameters
SF = 44000;
FeatureVector = zeros(size(waveforms,2),3);

%% Extract Features
% Normalize amplitudes
for i=1:size(waveforms,2)
waveforms(:,i) = waveforms(:,i)/max(waveforms(:,i)); 
end
%
for i=1:size(waveforms,2)
    [minval,minidx] = min(waveforms(:,i));
    [maxval,maxidx] = max(waveforms(minidx:end,i));
    FeatureVector(i,1) = maxidx/SF; % time from negative to positive peak
    maxidx = maxidx+minidx-1;
    FeatureVector(i,2) = abs(maxval)/abs(minval)/1000; % Ratio between positive to negative peak
    ReturnIdx = min(find(waveforms(maxidx:end,i)<=0));
    if isempty(ReturnIdx)
        [minval,ReturnIdx] = min(waveforms(maxidx:end,i));
    else
        FeatureVector(i,3) = ReturnIdx/SF; % time from positive peak to baseline
    end
end
%% Cluster wavefroms by features
[idx,C] = kmeans(FeatureVector,2,'Replicates',1000);

%% Plot Clustering Results
col=['r','g','b','m','c'];
for i=1:3
plot3(FeatureVector(idx==i,1),FeatureVector(idx==i,3),FeatureVector(idx==i,2),[col(i),'.'],'MarkerSize',12)
%plot(FeatureVector(idx==i,1)*1000,FeatureVector(idx==i,3)*1000,[col(i),'.'],'MarkerSize',12)
hold on
end
xlabel('Through to Peak Time')
ylabel('Peak to Baseline Time')
zlabel('Peak to Through Ratio')
%% Plot waveforms
t_waveforms = linspace(0,size(waveforms,1)/SF*10^3,size(waveforms,1)); 
% figure();
% for i=1:size(waveforms,2)
% plot(t_waveforms,waveforms(:,i),[col(idx(i)),'-'])
% hold on
% end
% Mean Waveforms
figure(); ax = axes();
for i=1:2
plot(t_waveforms,mean(waveforms(:,idx==i),2),[col(i),'-'],'LineWidth',2)
hold on
end
xlabel('Time[ms]')
ylabel('Scaled Amplitude')
xlim([0 3.5])
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
axes(ax)
end