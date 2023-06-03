%% Load from Figure
% Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
% SignalFiles=dir(fullfile(Dir, '*.fig'));
% pathname = Dir;
% prompt = {'1 = Baseline Activity subtraction | 0 = Show Baseline Activity'};
% dlgtitle = 'Input'; dims = [1 35]; definput = {'1'};
% answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
% 
% for i=1:length(SignalFiles)
%     openfig([pathname,'\',SignalFiles(i).name ]);
%     h=get(gca,'children');
%     intensity=get(h,'xdata');
%     Intensity{i} = intensity;
%     rate=get(h,'ydata');
%     if answer == '1'
%     %Rate{i} = (rate{end}-rate{1}(1))/max(rate{end}); % Subtract baseline activity from all responses and normalize to maximal response    
%     %Spon(i) = rate{1}(1)/max(rate{end});
%     Rate{i} = (rate{end}-rate{1}(1)); % Subtract baseline activity from all responses and normalize to maximal response    
%     Spon(i) = rate{1}(1);
%     
%     for k=1:length(Rate{i}) % Zero any negative values
%         if Rate{i}(k) < 0
%             Rate{i}(k) = 0;
%         end
%     end
%     Rate{i} = Rate{i}/max(Rate{i}); % Scale the fitted rate to the maximal response.
%     else
%     % Use this section for plots without baseline subtraction
%     Rate{i} = rate{end}/max(rate{end}); % Scale the fitted rate to the maximal response.
%     Spon(i) = rate{1}(1)/max(rate{end}); % Scale the baseline activity to the maximal response.
%     end
% pause()
% end

%% Load from Data file
clear all
close all
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = Dir;

prompt = {'1 = Baseline Activity subtraction | 0 = Show Baseline Activity'};
dlgtitle = 'Input'; dims = [1 35]; definput = {'1'};
answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
DurVec = strings;
for i=1:length(SignalFiles)
load([pathname,'\',SignalFiles(i).name ])
Intensity{i} = Data.NaturalIntensity;
    if answer == '1'
Rate{i} = Data.NaturalIntensityCount{1} - mean(Data.Spon{1});
for k=1:length(Rate{i}) % Zero any negative values
        if Rate{i}(k) < 0
            Rate{i}(k) = 0;
        end
end
    else
Rate{i} = Data.NaturalIntensityCount{1};
Spon{i} = Data.Spon{1};
    end
 %Rate{i} = Rate{i}/(Rate{i}(end)); % Scale the fitted rate to the maximal intensity.    
 Rate{i} = Rate{i}/max(Rate{i}); 
 DurVec = [DurVec;Data.StimDur];   
end
DurVec = DurVec(2:end);
%% Convert the Intensity Values to Log Scale
close all
for i=1:length(Intensity)
LogIntensity{i} = (Intensity{i});
end
%%
close all
for i=1:length(SignalFiles)
fig{i} = figure();
if answer == '1'
%plot(Intensity{i},Rate{i},'*')
%plot(LogIntensity{i},Rate{i},'*')
semilogx(LogIntensity{i},Rate{i},'*')
xlabel('Intensity Log [nW/mm^2]','FontSize',20)
ylabel('Scaled Spike Count','FontSize',20)
else
%plot(Intensity{i},Rate{i},'*',Intensity{i},ones(1,length(Intensity{i}))*Spon(i),'--')
semilogx(LogIntensity{i},Rate{i},'*')
xlabel('Intensity [nW/mm^2]','FontSize',15)
ylabel('Firing Rate [Spikes/s]','FontSize',15)
end
title(['Pulse Duration = ',DurVec{i},'ms']);
hold on
end
%% fit the data to a sigmoid  {\displaystyle f(x)={\frac {1}{1+e^{-x}}}}
for i=1:length(fig)
    figure(fig{i});
    % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}

%     fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*Intensity{i})))).^2); % Compute for Original Intensity Vals
    fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*LogIntensity{i})))).^2); % Compute for log Intensity Vals
    [x Eval(i)]=fminsearch(fun, [10 10 1]);
    %int_fit{i}=linspace(Intensity{i}(1),42.7,100); % Fit with minimal value presented
    int_fit{i}=linspace(0,max(LogIntensity{i}),1000); % Fit with 0 value sa minimum
    rate_fit{i} =(x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
    if answer =='0'
    SponNorm(i) = Spon(i);
    end
    hold on
    plot(int_fit{i},rate_fit{i})
    xlim([min(int_fit{i}) max(int_fit{i})])
    ylim([0 1.5])
    hold on
end
%% Build Intensities
IntVec = [];
for i=1:length(Intensity)
IntVec = [IntVec,Intensity{i}']; 
end
IntVec2 = unique(IntVec);
%% Select Relevant Units
RateVec = zeros(1,length(Rate));
for i = 1:length(Rate)
    if max(Rate{i})<2
    RateVec(i) = 1;
    end
end

Input = find(RateVec);
prompt = {'Select Relevant Units'};
dlgtitle = 'Input'; dims = [1 35]; definput = {num2str(Input)};
answer = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
RateNew = Rate(answer);
IntNew = Intensity(answer);
%% Build Response Per Intensity
RateMat = nan(length(RateNew),length(IntVec2));
for i = 1:length(RateNew)
for k = 1:length(RateNew{i})
   RateMat(i,find(IntVec2==IntNew{i}(k))) = RateNew{i}(k); 
end
end
NumDurs = size(unique(DurVec),1);
NewDurs = unique(DurVec);
NewDursVec = DurVec(answer);
for i = 1:NumDurs
    ind = find(NewDursVec == NewDurs{i});
    RateAvg{i} = mean(RateMat(ind,:),1,'omitnan');
    IntVec3{i} = IntVec2(find(~isnan(RateAvg{i})));
    RateAvgNew{i} = rmmissing(RateAvg{i});
end

%% Fit Avg Rate
    for i=1:NumDurs
    fun=@(X)sum((RateAvgNew{i}-(X(1)./(1+X(2)*exp(-X(3)*IntVec3{i})))).^2); % Compute for log Intensity Vals
    [X Eval(i)]=fminsearch(fun, [10 10 1]);
    %int_fit{i}=linspace(Intensity{i}(1),42.7,100); % Fit with minimal value presented
    int_fit_avg = linspace(0,max(IntVec3{i})+400,1000); % Fit with 0 value as minimum
    rate_fit_avg{i} = (X(1)./(1+X(2)*exp(-X(3)*int_fit_avg)));
    end
%% Find Stimulation Threshold
%Thresh = StimThresh(SponNorm);
%% Plot Scaled Fits of All Units
figure(); ax = axes();
color = ['r','b','g','h'];
count=1;
%     for k=1:length(IntVec2)
% %     ScaledRate(i,count) = rate_fit_avg(1,Idx(k));
%       Idx(k) = max(find(int_fit_avg<=IntVec2(k)));    
%       Int(count) = int_fit_avg(Idx(k));
%       RateForStd = rate_fit_avg(Idx(k));
%     count=count+1;
%     end
for i=1:length(NewDurs)
A(i) = str2num(NewDurs(i));
end
[B,C] = sort(A);


ScaledStd = std(RateMat,'omitnan')/sqrt(size(RateMat,2));
ScaledAvg = mean(RateMat,'omitnan');
for i=1:NumDurs
semilogx(int_fit_avg,rate_fit_avg{C(i)},color(i))
hold on
%errorbar(Int,rate_fit_avg(Idx),ScaledStd,'*','Color','b')
%hold on
%ScaledSpon = ones(1,length(ScaledAvg))*Thresh;
%plot(Int,ScaledSpon,'r--')
legends{i} = NewDurs(C(i));
hold on
end
legend([legends,'ms'],"Orientation","horizontal")
xlabel('Intensity [nW/mm^2]','FontSize',20)
ylabel('Scaled Spike Count Rate','FontSize',20)
xlim([-1 max(int_fit_avg)])
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
axes(ax)
%% Plot Raw Curves
figure();
color = ['r','b','g','h'];
for i=1:NumDurs
semilogx(IntVec3{i},RateAvgNew{i},color(i))
hold on
legends{i} = NewDurs(i);
end
legend([legends,'ms'],"Orientation","horizontal")
%% Dynamic Range Calc

IndMin = min(find(rate_fit_avg >= 0.1));
IndMax = min(find(rate_fit_avg >= 0.9));

ValDMin =  int_fit_avg(IndMin(i));
ValDMax =  int_fit_avg(IndMax(i));

DynamicRange =  20*log10(ValDMax)/ValDMin; %Calculate Dynamic Range in dB.
% bitDepth(i) = log(DynamicRange(i));
