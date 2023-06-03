
%% Load from Figures
clear all
close all
% Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
% SignalFiles=dir(fullfile(Dir, '*.fig'));
% pathname = Dir;
% Scale = [0.00042,0.02,0.27,0.3,0.9,1.01,1.45,1.47,1.93,2.09,2.38,2.8,3.03,3.15,3.6,4.05]*1000; % Intensity Vals in micro-Watts/mm^2
% RateMat{1} = zeros(length(SignalFiles),length(Scale)); RateMat{2} = zeros(length(SignalFiles),length(Scale));
% RateMatScale{1} = zeros(length(SignalFiles),length(Scale)); RateMatScale{2} = zeros(length(SignalFiles),length(Scale));
% UnitCount = zeros(1,2);
% 
% prompt = {'1 = Baseline Activity subtraction | 0 = Show Baseline Activity'};
% dlgtitle = 'Input';
% dims = [1 35];
% definput = {'1'};
% answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
% 
% MinInt = 0.5; MaxInt = 1;
% AllInt = [];     ResponseMat = NaN(length(SignalFiles),21);
% 
% for i=1:length(SignalFiles)
%     rate = {};
%     if ~isempty(regexp(SignalFiles(i).name,'10ms'))
%         Unit = 1;
%         RateMat{2}(i,:) = nan; RateMatScale{2}(i,:) = nan;
%         UnitCount(1) = UnitCount(1) + 1;
%     else
%         Unit = 2;
%         RateMat{1}(i,:) = nan; RateMatScale{1}(i,:) = nan;
%         UnitCount(2) = UnitCount(2) + 1;
%     end
% 
% 
% 
%     openfig([pathname,'\',SignalFiles(i).name ]);
%     %fig{i} = ans;
%     %title([SignalFiles(i).name]);
%     h=get(gca,'children');
%     intensity=get(h,'xdata');    
%     %Intensity{i} = intensity{1}*1000; %Convert Vals to micro-Watts/mm^2
%     
%     Intensity{i} = intensity{1};    
%     if min(Intensity{i})<MinInt
%         MinInt = min(Intensity{i});
%     end
%     if max(Intensity{i}) > MaxInt
%         MaxInt = max(Intensity{i});
%     end
%     AllInt = [AllInt,Intensity{i}];
%     UAllInt = unique(AllInt);
%     
% 
% 
%     rate=get(h,'ydata');
%        if answer == '1'
%     ResponseRate{i} = (rate{2}-rate{1}(1)); %Substracts Baseline Activity of Every Unit
%        end
% 
% 
%     %     if Intensity{i}(1) == 0
%     %         Intensity{i} = Intensity{i}(2:end);
%     %         ResponseRate{i} = ResponseRate{i}(2:end);
%     %     end
%     
%     % Zero any negative Vals after baseline subtraction
%     for  z = 1:length(ResponseRate{i})
% 
% 
%         if ResponseRate{i}(z) <0
%             ResponseRate{i}(z) = 0;
%         end
%     end
% 
%  
%         %Rate{i} = ResponseRate{i}/max(ResponseRate{i}); %  Normalize to Maximal Response.       
%         Rate{i} = ResponseRate{i}/ResponseRate{i}(end);
%         
%         %Rate{i} = ResponseRate{i};
% %     else
% %         Rate{i} = rate{2}/max(rate{2});
% %         Spon(i) = rate{1}(1)/max(rate{2});
%     
% 
% 
% 
%     UnitIdx(i) = Unit;
%     pause();
% end
% 
% for i = 1:length(Intensity)
%     for k=1:length(Intensity{i})
%         A = find(UAllInt == Intensity{i}(k));
%         ResponseMat(i,A) = Rate{i}(k);
%     end
% end
%% Load from Data files
clear all
close all
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = Dir;

prompt = {'1 = Baseline Activity subtraction | 0 = Show Baseline Activity'};
dlgtitle = 'Input'; dims = [1 35]; definput = {'1'};
answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
DurVec = strings; IntVec = [];
for i=1:length(SignalFiles)
load([pathname,'\',SignalFiles(i).name ])
Intensity{i} = Data.ProstheticIntensity;
IntVec = [IntVec;Data.ProstheticIntensity];
if answer == '1'
Rate{i} = Data.ProstheticIntensityCount{1} - (cell2mat(Data.Spon)*4)'; % Remove Baseline of each trial
%Rate{i} = Data.ProstheticIntensityCount{1} - mean(cell2mat(Data.Spon))*4; % Remove Mean Baseline from all trials.
for k=1:length(Rate{i}) % Zero any negative values
        if Rate{i}(k) < 0
            Rate{i}(k) = 0;
        end
end
    else
Rate{i} = Data.ProstheticIntensityCount{1};
Spon{i} = Data.Spon{1};
    end
 Rate{i} = Rate{i}/(Rate{i}(end)); % Scale the fitted rate to the maximal intensity.    
 DurVec = [DurVec;Data.StimDur];   
end
DurVec = DurVec(2:end);

%% Build Response Mat 
AllIntensity = unique(IntVec); RateMat = nan(length(SignalFiles),length(AllIntensity));
for r=1:length(SignalFiles)
Ind = find(ismember(AllIntensity,Intensity{r}));
RateMat(r,Ind) = Rate{r};
end
 %% Convert the Intensity Values to Log Scale
% for i=1:length(Intensity)
%     LogIntensity{i} = log10(Intensity{i});
% end
%% Plot Dots
close all
for i=1:length(SignalFiles)
    fig{i} = figure();
    if answer == '1'
        semilogx(Intensity{i}*1000,Rate{i},'*')
        %xlim([0, max(Intensity{i})+1]);
    else
        plot(LogIntensity{i},Rate{i},'*',Intensity{i},ones(1,length(Intensity{i}))*Spon(i),'--')
    end
end
%% Fit Single Units
% RateAvg{1} = mean([Rate{1};Rate{2}],1); % 4ms
% RateAvg{2} = mean([Rate{3};Rate{4}],1); % 10ms
% RateAvg{3} = mean([Rate{1};Rate{2};Rate{3}(:,1:end-2);Rate{4}(:,1:end-2)]); % all
%IntensityAvg{1} = Intensity{1}; IntensityAvg{2} = Intensity{3}; IntensityAvg{3} = Intensity{1};
for i=1:length(fig)
    figure(fig{i});
    % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}
    fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*Intensity{i})))).^2);
    [x Eval(i)]=fminsearch(fun, [10 10 100]);
    %options = optimset('MaxFunEvals',10000,'MaxIter',10000);
    %[x Eval(i)]=fminsearch(fun, [max(Rate{i}) 1 1],options);
    int_fit{i}=linspace(min(Intensity{i}),max(Intensity{i}),1000);
    rate_fit{i} =(x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
    %rate_fit_norm{i} = (rate_fit{i})/max(rate_fit{i}); % Scale the fitted rate to the maximal response.
    %SponNorm(i) = Spon(i)/max(rate_fit{i}); % Scale the baseline activity to the maximal response.
    rate_fit_norm{i} = (rate_fit{i});
    %SponNorm(i) = Spon(i);
    hold on
    plot(int_fit{i}*1000,rate_fit{i})
    set(gca,'xscale','log')
    %xlim([1 max(int_fit{i}*1000)])
    xlabel('Intensity [\muW/mm^2]')
    ylabel('Scaled Firing Rate')
end
%% Single Unit Dynamic Range
% SelectUnits = [1,2,3,5,6];
% 
% for i = 1:length(SelectUnits)
% IndMin(i) = min(find(rate_fit{SelectUnits(i)} >= 0.1));
% IndMax(i) = min(find(rate_fit{SelectUnits(i)} >= 0.9));
% 
% ValDMin(i) =  int_fit{i}(IndMin(i))*1000;
% ValDMax(i) =  int_fit{i}(IndMax(i))*1000;
% 
% DynamicRange(i) =  20*log10(ValDMax(i)/ValDMin(i)); %Calculate Dynamic Range in dB.
% end
%%
RateVec = zeros(1,length(Rate));
for i = 1:length(Rate)
    if max(Rate{i})<3 % Select Only Units with small variations in the response
        RateVec(i) = 1;
    end
end

Input = find(RateVec);

prompt = {'Select Units to include in further analysis'};
dlgtitle = 'Input';
dims = [1 35];
%definput = {'[1,2,3,4,5,6,7,8,9,10,11,12]'};
definput = {num2str(Input)};
answer1 = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
StimDur = nan(1,length(SignalFiles));
for k=1:length(SignalFiles)
    %DurIdx = strfind(SignalFiles(answer1(k)).name,'ms');
    %StimDur(k) = str2num(SignalFiles(answer1(k)).name(DurIdx-2:DurIdx-1));
    if ismember(k,answer1)
    StimDur(k) = DurVec(k);
    end
end

Durs = rmmissing(unique(StimDur)); % Count num of different durations
CountMat = cell(2,1);
for k=1:length(Durs)
%     for i = 1:length(answer1)
%         if StimDur(i) == Durs(k)
%             CountMat{k} = [CountMat{k};rate_fit{answer1(k)}];
%         end
%     end
DurIdx{k} =  find(ismember(StimDur,Durs(k)));
end

% Average across Stim durs in all units and extract relevant Intensities
prompt = {'Insert Idxs for 4ms:','Insert Idxs for 10ms:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {num2str(DurIdx{1});num2str(DurIdx{2})};
answer2 = inputdlg(prompt,dlgtitle,dims,definput);

% for i=1:length(answer2) 
% AvgCount(i,:) = mean(ResponseMat(min(str2num(answer2{i})):max(str2num(answer2{i})),:),'omitnan');
% StdCount(i,:) = std(ResponseMat(min(str2num(answer2{i})):max(str2num(answer2{i})),:),'omitnan');
% AvgInt{i} = UAllInt(find( ~isnan(AvgCount(i,:)))); 
% end
for i = 1:length(Durs)
AvgCount(i,:) = mean(RateMat(find(StimDur == Durs(i)),:),'omitnan');
StdCount(i,:) = std(RateMat(find(StimDur == Durs(i)),:),'omitnan');
AvgInt{i} = AllIntensity(find( ~isnan(AvgCount(i,:))))';
end

%% Fit Average Results
for k = 1:size(AvgCount,1)
    fun=@(x)sum((rmmissing(AvgCount(k,:))-(x(1)./(1+x(2)*exp(-x(3)*AvgInt{k})))).^2);
    [x Eval(i)]=fminsearch(fun, [1 1 100]);
    AVGint_fit{k}=linspace(min(AvgInt{k}),max(AvgInt{k}),1000);
    AVGrate_fit{k} =(x(1)./(1+x(2)*exp(-x(3)*AVGint_fit{k})));
end

%% Plot Average Results
col = ['r','b'];
IdxError = [9,14,17];
figure();
ax = axes();
for k = 1:length(Durs)
h(k) = semilogx(AVGint_fit{k}*1000,AVGrate_fit{k},'Color',col(k));
%legend('AutoUpdate','off')
hold on
[a,b] = intersect(round(AVGint_fit{k},2),round(AllIntensity(IdxError),2),'stable');
%errorbar(AVGint_fit{k}(b)*1000,AVGrate_fit{k}(b),StdCount(k,IdxError)/sqrt(size(RateMat,2)),'.','Color',col(k)) % Plot SEM 
%errorbar(AVGint_fit{k}(b)*1000,AVGrate_fit{k}(b),StdCount(k,IdxError),'.','Color',col(k)) % Plot STD
%plot(AvgInt{k}*1000,rmmissing(AvgCount(k,:)),'*','Color',col(k));
hold on
end
legend([h(1) h(2)],['4ms ';'10ms'])
%xlim([10 max(AVGint_fit{k}*1000)])
%ylim([0 1.1])
xlabel('Intensity [\muW/mm^2]')
ylabel('Scaled Spike Count')
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
%ax.XLim = [0,3500];
axes(ax)

