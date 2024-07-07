clear all
close all
% Define path for data
Dir=uigetdir('*.fig','Select a Folder to Load Data Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = [Dir,'\'];

% Define variables
IntNatural = []; IntProsthetic = [];
LatencyNatural = []; LatencyProsthetic = [];
waveforms = []; Idxwaves = []; NumSpikes = [];
LatencyforSpikes = {};
% Iterate over data files
countProsthetic = 1; countnatural = 1;
for i= 1:length(SignalFiles)
    load([pathname,SignalFiles(i).name])
    if isfield(Data,'NaturalLatency')
        LatencyNatural = [LatencyNatural; cell2mat(Data.NaturalLatency)];
        for k = 1:length(Data.NaturalIntensity)
            if ~isnan(Data.NaturalLatency{1}(k))
                IntNatural = [IntNatural,Data.NaturalIntensity(k)];
            end
        end
        countnatural = countnatural+1;
        waveforms = [waveforms,Data.SortedSpikes{1}{Data.Clusters}];
        Idxwaves = [Idxwaves,ones(1,size(Data.SortedSpikes{1}{Data.Clusters},2))];
        NumSpikes = [NumSpikes,size(Data.SortedSpikes{1}{Data.Clusters},2)];
        A = rmmissing(cell2mat(Data.NaturalLatency));
        if length(A)>1
            LatencyforSpikes{i} = A(end-1:end);
        else
            LatencyforSpikes{i} = A;
        end
    else
        LatencyProsthetic = [LatencyProsthetic; cell2mat(Data.ProstheticLatency)];
        for k = 1:length(Data.ProstheticIntensity)
            if ~isnan(Data.ProstheticLatency{1}(k))
                IntProsthetic = [IntProsthetic,Data.ProstheticIntensity(k)];
            end
        end
        countProsthetic = countProsthetic+1;
        waveforms = [waveforms,Data.SortedSpikes{1}{Data.Clusters}];
        Idxwaves = [Idxwaves,ones(1,size(Data.SortedSpikes{1}{Data.Clusters},2))*2];
        NumSpikes = [NumSpikes,size(Data.SortedSpikes{1}{Data.Clusters},2)];
        A = rmmissing(cell2mat(Data.ProstheticLatency));
        if length(A)>1
            LatencyforSpikes{i} = A(end-1:end);
        else
            LatencyforSpikes{i} = A;
        end
    end
    close all
end

%% Calc Mean Latency Per Intensity
LatencyNatural = rmmissing(LatencyNatural); LatencyProsthetic = rmmissing(LatencyProsthetic);
AllIntP = unique(IntProsthetic); AllIntN = unique(IntNatural);

LatMatN = nan([length(AllIntN) length(LatencyNatural)]);countN = 1;
for i = 1:length(AllIntN)
    LatVecN(i) = mean(LatencyNatural(find(AllIntN(i) == IntNatural)));


    for k = 1:length(LatencyNatural)
        if IntNatural(k) == AllIntN(i)
            LatMatN(i,countN) = LatencyNatural(k);
            countN = countN +1;
        else
            LatMatN(i,countN) = nan;
        end
    end
end
LatN = mean(LatMatN,2,"omitnan");

LatMatP = nan([length(AllIntP) length(LatencyProsthetic)]);countP = 1;
for i = 1:length(AllIntP)
    for k = 1:length(LatencyProsthetic)
        if IntProsthetic(k) == AllIntP(i)
            LatMatP(i,countP) = LatencyProsthetic(k);
            countP = countP +1;
        else
            LatMatP(i,countP) = nan;
        end
    end
end
LatP = mean(LatMatP,2,"omitnan");

%% Calc Total Latency
AvgLNat = mean(LatencyNatural,'omitnan'); AvgLPros = mean(LatencyProsthetic,'omitnan');
StdLNat = std(LatencyNatural,'omitnan'); StdLPros = std(LatencyProsthetic,'omitnan');

%% Plot
FN = fit(AllIntN',LatN,'poly1');
Nfig = figure();axN = axes();
plot(FN,AllIntN,LatN)
ylabel('Latency [ms]')
xlabel ('Intensity [nW/mm^2]')
ylim([0 120])
axN.PlotBoxAspectRatio = [1,1,1]; axN.FontSize = 20;
axN.Box = 'off'; axN.Color = "none"; legend('off')
axes(axN)
%set(axN,'XScale','log')
%xlabel(log10(AllIntN))
Pidx = find(~isnan(LatP));
FP = fit(AllIntP(Pidx)',LatP(Pidx),'poly1');
Pfig = figure();axP = axes();
plot(FP,AllIntP,LatP)
ylabel('Latency [ms]')
xlabel ('Intensity [mW/m^2]')
%ylim([0 120])
axP.PlotBoxAspectRatio = [1,1,1]; axP.FontSize = 20;
axP.Box = 'off'; axP.Color = "none"; legend('off')
axes(axP)
%set(axP,'XScale','log')

%%
figure(); ax = axes();
Names = categorical({'Natural (N =24)','Prosthetic (N=12)'});
bar(Names,[AvgLNat,AvgLPros],0.4)
hold on
errorbar(Names,[AvgLNat,AvgLPros],[StdLNat,StdLPros],'.','LineWidth',2)
ylabel('Latency [ms]')
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
axes(ax)
%% Cluster Waveforms
[ClustIdx,color] = WaveClust(waveforms,Idxwaves);

%% plot Latecy distrinbutions
count = 1;
while strfind(SignalFiles(count).name,'Natural') == 1
    count=count+1;
end
NumNUnits = count-1; NumPUnits = size(SignalFiles,1) - NumNUnits;

VecN = zeros(2,150); VecP = zeros(2,150); spikecount = 1;
for i = 1:length(VecN)
    for k=1:NumNUnits
        X = unique(LatencyforSpikes{k});
        for c =1:length(X)
            VecN(1,X) = VecN(1,X)+length(find(ClustIdx(spikecount:spikecount+NumSpikes(k)-1)==1));
            VecN(2,X) = VecN(2,X)+length(find(ClustIdx(spikecount:spikecount+NumSpikes(k)-1)==2));
            spikecount = sum(NumSpikes(1:k));
        end
    end
end
for i = 1:length(VecP)
    for k=NumNUnits+1:NumNUnits+NumPUnits
        X = unique(LatencyforSpikes{k});
        for c =1:length(X)
            VecP(1,X) = VecP(1,X)+length(find(ClustIdx(spikecount:spikecount+NumSpikes(k)-1)==1));
            VecP(2,X) = VecP(2,X)+length(find(ClustIdx(spikecount:spikecount+NumSpikes(k)-1)==2));
            spikecount = sum(NumSpikes(NumNUnits+1:k));
        end
    end
end
x = find(VecN(1,:));
VecN1 = VecN(1,x)/2; VecN2 = VecN(2,x)/2;
z = find(VecP(1,:));
VecP1 = VecP(1,z)/2; VecP2 = VecP(2,z)/2;
%%
figure();
subplot(1,2,1)
bar(x,VecN1,'FaceColor','red')
hold on
bar(x,VecN2,'FaceColor','green')
xlabel('Latency [ms]'); ylabel('Spike Count'); title('Natural');

subplot(1,2,2)
bar(z,VecP1,'FaceColor','red')
Fontsize = 20;
hold on
bar(z,VecP2,'FaceColor','green')
xlabel('Latency [ms]'); ylabel('Spike Count'); title('Prosthetic');

%%
figure();
subplot(1,2,1)
histogram(LatMatN(end-1:end,:)',14)
xlabel('Latency [ms]'); ylabel('Count');
xlim([0 130]); ylim([0 10]);
title('Natural')
hold on
subplot(1,2,2)
histogram(LatMatP(end-1:end,:)',4)
xlabel('Latency [ms]'); ylabel('Count');
xlim([0 130]); ylim([0 30]);
title('Prosthetic')
