close all
clear all
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.fig'));
pathname = Dir;

prompt = {'1 = Baseline Activity subtraction | 0 = Show Baseline Activity'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1'};
answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));

for i=1:length(SignalFiles)
    openfig([pathname,'\',SignalFiles(i).name ]);
    h=get(gca,'children');
    intensity=get(h,'xdata');
    Intensity{i} = intensity{1};
    rate=get(h,'ydata');
    if answer == '1'
    %Rate{i} = (rate{end}-rate{1}(1))/max(rate{end}); % Subtract baseline activity from all responses and normalize to maximal response    
    %Spon(i) = rate{1}(1)/max(rate{end});
    Rate{i} = (rate{end}-rate{1}(1)); % Subtract baseline activity from all responses and normalize to maximal response    
    Spon(i) = rate{1}(1);
    
    for k=1:length(Rate{i}) % Zero any negative values
        if Rate{i}(k) < 0
            Rate{i}(k) = 0;
        end
    end
    Rate{i} = Rate{i}/max(Rate{i}); % Scale the fitted rate to the maximal response.
    else
    % Use this section for plots without baseline subtraction
    Rate{i} = rate{end}/max(rate{end}); % Scale the fitted rate to the maximal response.
    Spon(i) = rate{1}(1)/max(rate{end}); % Scale the baseline activity to the maximal response.
    end
pause()
end

%% Convert the Intensity Values to Log Scale
close all
for i=1:length(Intensity)
LogIntensity{i} = (Intensity{i});
end
%%
for i=1:length(SignalFiles)
fig{i} = figure();
if answer == '1'
%plot(Intensity{i},Rate{i},'*')
plot(LogIntensity{i},Rate{i},'*')
semilogx(LogIntensity{i},Rate{i},'*')
xlabel('Intensity Log [cd/m^2]','FontSize',15)
ylabel('Scaled Firing Rate','FontSize',15)
else
%plot(Intensity{i},Rate{i},'*',Intensity{i},ones(1,length(Intensity{i}))*Spon(i),'--')
semilogx(LogIntensity{i},Rate{i},'*',LogIntensity{i},ones(1,length(LogIntensity{i}))*Spon(i),'--')
xlabel('Intensity [cd/m^2]','FontSize',15)
ylabel('Firing Rate [Spikes/s]','FontSize',15)
end

end

%% fit the data to a sigmoid  {\displaystyle f(x)={\frac {1}{1+e^{-x}}}}

for i=1:length(fig)
    figure(fig{i});
    % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}

%     fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*Intensity{i})))).^2); % Compute for Original Intensity Vals
    fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*LogIntensity{i})))).^2); % Compute for log Intensity Vals
    [x Eval(i)]=fminsearch(fun, [1 1 1]);
    %int_fit{i}=linspace(Intensity{i}(1),42.7,100); % Fit with minimal value presented
    int_fit{i}=linspace(0,max(LogIntensity{i})+40,100); % Fit with 0 value sa minimum
    rate_fit{i} =(x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
    rate_fit_norm{i} = (rate_fit{i}/max(rate_fit{i})); % Normalize the Vals by maximal response of each unit 
    if answer =='0'
    SponNorm(i) = Spon(i);
    end
    hold on
    plot(int_fit{i},rate_fit{i})
end
%% Build Intensities
IntVec = [];
for i=1:length(Intensity)
IntVec = [IntVec,Intensity{i}]; 
end
IntVec2 = unique(IntVec);
%% Build Response Per Intensity
RateMat = nan(length(Rate),length(IntVec2));
for i = 1:length(Rate)
for k = 1:length(Rate{i})
   RateMat(i,find(IntVec2==Intensity{i}(k))) = Rate{i}(k); 
end
end
RateAvg = mean(RateMat,1,'omitnan');
%% Fit Avg Rate
    fun=@(x)sum((RateAvg-(x(1)./(1+x(2)*exp(-x(3)*IntVec2)))).^2); % Compute for log Intensity Vals
    [x Eval(i)]=fminsearch(fun, [1 1 1]);
    %int_fit{i}=linspace(Intensity{i}(1),42.7,100); % Fit with minimal value presented
    int_fit_avg=linspace(0,max(IntVec2)+40,100); % Fit with 0 value sa minimum
    rate_fit_avg =(x(1)./(1+x(2)*exp(-x(3)*int_fit_avg)));

%% Find Stimulation Threshold
%Thresh = StimThresh(SponNorm);
%% Plot Scaled Fits of All Units


figure(); ax = axes();
count=1;
    for k=1:length(IntVec2)
%     ScaledRate(i,count) = rate_fit_avg(1,Idx(k));
      Idx(k) = max(find(int_fit_avg<=IntVec2(k)));    
      Int(count) = int_fit_avg(Idx(k));
      RateForStd = rate_fit_avg(Idx(k));
    count=count+1;
    end

ScaledStd = std(RateMat,'omitnan')/sqrt(size(RateMat,2));
ScaledAvg = mean(RateMat,'omitnan');
semilogx(int_fit_avg,rate_fit_avg,'r')
hold on
errorbar(Int,rate_fit_avg(Idx),ScaledStd,'*','Color','b')
hold on
%ScaledSpon = ones(1,length(ScaledAvg))*Thresh;
%plot(Int,ScaledSpon,'r--')
ylim([0 1.1]);
xlabel('Intensity Log [cd/m^2]','FontSize',20)
ylabel('Scaled Spike Count Rate','FontSize',20)
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
axes(ax)
%% Dynamic Range Calc

IndMin = min(find(rate_fit_avg >= 0.1));
IndMax = min(find(rate_fit_avg >= 0.9));

ValDMin =  int_fit_avg(IndMin(i));
ValDMax =  int_fit_avg(IndMax(i));

DynamicRange =  20*log10(ValDMax)/ValDMin; %Calculate Dynamic Range in dB.
% bitDepth(i) = log(DynamicRange(i));
