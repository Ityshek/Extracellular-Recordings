clear all
close all
% Define path for data
Dir=uigetdir('*.fig','Select a Folder to Load Data Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = [Dir,'\'];

% Define variables
IntNatural = []; IntProsthetic = [];
LatencyNatural = []; LatencyProsthetic = [];
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
    else
        LatencyProsthetic = [LatencyProsthetic; cell2mat(Data.ProstheticLatency)];
        for k = 1:length(Data.ProstheticIntensity)
            if ~isnan(Data.ProstheticLatency{1}(k))
                IntProsthetic = [IntProsthetic,Data.ProstheticIntensity(k)];
            end
        end
        countProsthetic = countProsthetic+1;
    end
    close all
end

%% Calc Mean Latency Per Intensity
LatencyNatural = rmmissing(LatencyNatural); LatencyProsthetic = rmmissing(LatencyProsthetic);
AllIntP = unique(IntProsthetic); AllIntN = unique(IntNatural);

LatMatN = nan([length(AllIntN) length(LatencyNatural)]);countN = 1;
for i = 1:length(AllIntN)
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
xlabel ('Intensity [cd/m^2]')
ylim([0 120])
axN.PlotBoxAspectRatio = [1,1,1]; axN.FontSize = 20;
axN.Box = 'off'; axN.Color = "none"; legend('off')
axes(axN)
%set(axN,'XScale','log')
%xlabel(log10(AllIntN))

FP = fit(AllIntP',LatP,'poly1');
Pfig = figure();axP = axes();
plot(FP,AllIntP,LatP)
ylabel('Latency [ms]')
xlabel ('Intensity [mW/m^2]')
ylim([0 120])
axP.PlotBoxAspectRatio = [1,1,1]; axP.FontSize = 20;
axP.Box = 'off'; axP.Color = "none"; legend('off')
axes(axP)
%set(axP,'XScale','log')



%%
figure(); ax = axes();
Names = categorical({'Natural (N = 4)','Prosthetic (N = 7)'});
bar(Names,[AvgLNat,AvgLPros],0.4)
hold on
errorbar(Names,[AvgLNat,AvgLPros],[StdLNat,StdLPros],'.','LineWidth',2)
ylabel('Latency [ms]')
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
axes(ax)