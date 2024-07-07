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

Dates = strings; Channels = strings;
Locs = strings; DurVec = strings;

for i=1:length(SignalFiles)
    load([pathname,'\',SignalFiles(i).name ])
    Intensity{i} = ceil(Data.NaturalIntensity./10) * 10;
    if answer == '1'
        rate{i} = Data.NaturalIntensityCount{1} - mean(cell2mat(Data.Spon));
        for k=1:length(rate{i}) % Zero any negative values
            if rate{i}(k) < 0
                rate{i}(k) = 0;
            end
        end
    else
        rate{i} = Data.NaturalIntensityCount{1};
        Spon{i} = Data.Spon{1};
    end
    % Sort Trials by Unit
    [a,b] = regexp(SignalFiles(i).name,'Hz_[0123456789.]*');
    SortVar(i).Date = SignalFiles(i).name(a+3:b-1);
    [c,d] = regexp(SignalFiles(i).name,'Loc[0123456789]*');
    SortVar(i).Loc = SignalFiles(i).name(d);
    [e,f] = regexp(SignalFiles(i).name,'Ch[0123456789]*');
    SortVar(i).Channel = SignalFiles(i).name(e+2:f);
    SortVar(i).Duration = Data.StimDur;
    Dates = [Dates,string(SortVar(i).Date)];
    Channels = [Channels,string(SortVar(i).Channel)];
    Locs = [Locs,string(SortVar(i).Loc)];
    DurVec = [DurVec;Data.StimDur];
    clear a b c d e f
end
DurVec = DurVec(2:end); Dates = Dates(2:end);
Locs = Locs(2:end); Channels = Channels(2:end);

%
%Rate{i} = Rate{i}/(Rate{i}(end)); % Scale the fitted rate to the maximal intensity.
%Rate{i} = Rate{i}/max(Rate{i});

%% Choose Which Trials to Normalize
IdxPostNorm = []; % Contains idxs of the trials that can be normalized. 
for s=1:length(SortVar)
    DateMatch = find(matches(Dates,Dates(s))); % Find All Dates same as Current Date
    LocMatch = find(matches(Locs,Locs(s))); % Find All Locs same as Current Loc
    ChannelMatch = find(matches(Channels,Channels(s))); % Find All Channels same as Current Channel
    Match1 = intersect(DateMatch,LocMatch);
    Match2 = intersect(Match1,ChannelMatch);
    RelevantDurs = [char(DurVec(Match2))];
    [a,b] = max(str2num(RelevantDurs));
    if a == 10
        IdxPostNorm(end+1) = s;
        for v = 1:length(Match2)
            %Rate{Match2(v)} = rate{Match2(v)}/rate{Match2(b)}(end);
            Rate{Match2(v)} = rate{Match2(v)}/max(rate{Match2(b)});
        end
    else
        for v = 1:length(Match2)
            Rate{Match2(v)} = NaN;
        end
    end
end
%% Convert the Intensity Values to Log Scale
close all
for i=1:length(IdxPostNorm)
    LogIntensity{i} = (Intensity{IdxPostNorm(i)});
end
%% Plot Data Points
close all
for i=1:length(IdxPostNorm)
    fig{i} = figure();
    if answer == '1'
        %plot(Intensity{i},Rate{i},'*')
        %plot(LogIntensity{i},Rate{i},'*')
        semilogx(LogIntensity{i},Rate{IdxPostNorm(i)},'*')
        xlabel('Intensity Log [nW/mm^2]','FontSize',20)
        ylabel('Scaled Spike Count','FontSize',20)
    else
        %plot(Intensity{i},Rate{i},'*',Intensity{i},ones(1,length(Intensity{i}))*Spon(i),'--')
        semilogx(LogIntensity{i},Rate{IdxPostNorm(i)},'*')
        xlabel('Intensity [nW/mm^2]','FontSize',15)
        ylabel('Firing Rate [Spikes/s]','FontSize',15)
    end
    title(['Pulse Duration = ',DurVec{i},'ms']);
    hold on
end
%% fit the data to a sigmoid  {\displaystyle f(x)={\frac {1}{1+e^{-x}}}}
for i=1:length(IdxPostNorm)
    figure(fig{i});
    % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}

    %     fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*Intensity{i})))).^2); % Compute for Original Intensity Vals
    %fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*LogIntensity{i})))).^2); % Compute for log Intensity Vals
    fun=@(x)sum((Rate{IdxPostNorm(i)}-(x(1)./(1+x(2)*exp(-x(3)*LogIntensity{i})))).^2); % Compute for log Intensity Vals
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
    %ylim([0 1.5])
    hold on
end

%% Select Relevant Units
RateVec = zeros(1,length(IdxPostNorm));
for i = 1:length(IdxPostNorm)
    if max(Rate{IdxPostNorm(i)})<3
        RateVec(i) = 1;
    end
%     if max(rate_fit{i}) < 1.5
%         RateVec(i) = 1;
%     end
end
%RateVec = [1 2 3 4 5 6 7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
Input = find(RateVec);
prompt = {'Select Relevant Units'};
dlgtitle = 'Input'; dims = [1 35]; definput = {num2str(Input)};
answer = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
RateNew = Rate(IdxPostNorm(answer));
IntNew = Intensity(IdxPostNorm(answer));
%% Build Intensities
IntVec = [];
for i=1:length(Intensity)
    IntVec = [IntVec,Intensity{i}'];
end
IntVec2 = unique(IntVec);
%% Build Response Per Intensity
RateMat = nan(length(RateNew),length(IntVec2));
for i = 1:length(RateNew)
    for k = 1:length(RateNew{i})
        RateMat(i,find(IntVec2==IntNew{i}(k))) = RateNew{i}(k);
        %RateMat(i,find(IntVec2==Intensity{IdxPostNorm(i)}(k))) = Rate{IdxPostNorm(i)}(k);
    end
end
NumDurs = size(unique(DurVec),1);
NewDurs = unique(DurVec);
NewDursVec = DurVec(IdxPostNorm);
for i = 1:NumDurs
    ind{i} = find(NewDursVec == NewDurs{i});
    RateAvg{i} = mean(RateMat(ind{i},:),1,'omitnan');
    RateSEM{i} = std(RateMat(ind{i},:),1,'omitnan')/sqrt(length(RateMat(ind{i})));
    IntVec3{i} = IntVec2(find(~isnan(RateAvg{i})));
    RateAvgNew{i} = rmmissing(RateAvg{i});
end

%% Fit Avg Rate
% Fit as it should be..
for i=1:NumDurs
    fun=@(X)sum((RateAvgNew{i}-(X(1)./(1+X(2)*exp(-X(3)*IntVec3{1})))).^2); % Compute for log Intensity Vals
    %options = optimset('MaxFunEvals',100000000,'MaxIter',100000000,'TolFun',2e-8);    
    [X Zval(i)]=fminsearch(fun, [1 1 1]);
    %int_fit{i}=linspace(Intensity{i}(1),42.7,100); % Fit with minimal value presented
    int_fit_avg = linspace(0,10000,10000); % Fit with 0 value as minimum
    rate_fit_avg{i} = (X(1)./(1+X(2)*exp(-X(3)*int_fit_avg)));
end
%Find Fits forced to 0 and 1
% for i=1:NumDurs
%     fun=@(X)sum((RateAvgNew{i}-(1./(1+1*exp(-X(1)*IntVec3{i})))).^2); % Compute for log Intensity Vals
%     options = optimset('MaxFunEvals',100000000,'MaxIter',100000000,'TolFun',2e-8);    
%     [X Zval(i)]=fminsearch(fun, [1],options);
%     %int_fit{i}=linspace(Intensity{i}(1),42.7,100); % Fit with minimal value presented
%     int_fit_avg = linspace(0,max(IntVec3{i}),10000); % Fit with 0 value as minimum
%     rate_fit_avg{i} = (1./(1+1*exp(-X(1)*int_fit_avg)));
% end   
%%
save('C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Results\Data\Natural\Post-Fit\Data_10msMaxIntensityScale2.mat')
%% Plot Scaled Fits of All Units
figure(); ax = axes();
color = ['r','b','g']; Line = ['-',':','--'];
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
    h{i} = semilogx(int_fit_avg,rate_fit_avg{C(i)},'Color',color(C(i)),'LineWidth',2,LineStyle=Line(C(i)));
    hold on
    SEMY = rate_fit_avg{C(i)}(int32(IntVec3{1}(6:end)));
    errorbar(IntVec3{1}(6:end),SEMY,rmmissing(RateSEM{C(i)}),'LineStyle','none','Color',color(C(i)))
    %hold on
    %ScaledSpon = ones(1,length(ScaledAvg))*Thresh;
    %plot(Int,ScaledSpon,'r--')
    legends{i} = NewDurs(C(i));
    hold on
end

% Plot Dots on Average Fits
% for i=1:length(ind) % loop over all durs
%     a = ind{i}; % Convert relevant unit idxs in current stim dur to nums
%     for k = 1:length(a)
%         hold on
%         if i==1
%             semilogx(IntNew{a(k)},rmmissing(RateMat(a(k),:)),'*r')
%         elseif i==2
%             semilogx(IntNew{a(k)},rmmissing(RateMat(a(k),:)),'*b')
%         else
%             semilogx(IntNew{a(k)},rmmissing(RateMat(a(k),:)),'*g')
%         end
%     end
% end
legend([h{1},h{2},h{3}],{'2ms';'4ms';'10ms'},"Orientation","vertical","Location","southeast","Box","off")
xlabel('Intensity [nW/mm^2]','FontSize',20)
ylabel('Scaled Spike Count','FontSize',20)
title('Natural Vision')
xlim([1 max(int_fit_avg)]); ylim([0 1.5])
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
axes(ax);
%% Plot Raw Curves
figure();
color = ['r','b','g','h'];
for i=1:NumDurs
    semilogx(IntVec3{i},RateAvgNew{i},color(i))
    hold on
    legends{i} = NewDurs(i);
end
legend([legends,'ms'],"Orientation","horizontal")
%% Find Stimulation Threshold
%Thresh = StimThresh(SponNorm);
%% Dynamic Range Calc

IndMin = min(find(rate_fit_avg >= 0.1));
IndMax = min(find(rate_fit_avg >= 0.9));

ValDMin =  int_fit_avg(IndMin(i));
ValDMax =  int_fit_avg(IndMax(i));

DynamicRange =  20*log10(ValDMax)/ValDMin; %Calculate Dynamic Range in dB.
% bitDepth(i) = log(DynamicRange(i));
