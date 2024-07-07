
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
% for i=1:length(SignalFiles)
% load([pathname,'\',SignalFiles(i).name ])
% Intensity{i} = Data.ProstheticIntensity;
% IntVec = [IntVec;Data.ProstheticIntensity];
% if answer == '1'
% Rate{i} = Data.ProstheticIntensityCount{1} - (cell2mat(Data.Spon)*4)'; % Remove Baseline of each trial
% %Rate{i} = Data.ProstheticIntensityCount{1} - mean(cell2mat(Data.Spon))*4; % Remove Mean Baseline from all trials.
% for k=1:length(Rate{i}) % Zero any negative values
%         if Rate{i}(k) < 0
%             Rate{i}(k) = 0;
%         end
% end
%     else
% %Rate{i} = Data.ProstheticIntensityCount{1};
% Spon{i} = Data.Spon{1};
%     end
%  %Rate{i} = Rate{i}/(Rate{i}(end)); % Scale the fitted rate to the maximal intensity.    
%  DurVec = [DurVec;Data.StimDur];   
% end

Dates = strings; Channels = strings;
Locs = strings; DurVec = strings;

for i=1:length(SignalFiles)
    load([pathname,'\',SignalFiles(i).name ])
    intensity{i} = Data.ProstheticIntensity; IntVec = [IntVec;Data.ProstheticIntensity];
    if answer == '1'
        rate{i} = Data.ProstheticIntensityCount{1} - mean(cell2mat(Data.Spon));
        for k=1:length(rate{i}) % Zero any negative values
            if rate{i}(k) < 0
                rate{i}(k) = 0;
            end
        end
    else
        rate{i} = Data.ProstheticIntensityCount{1};
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

close all
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
            Rate{Match2(v)} = rate{Match2(v)}/rate{Match2(b)}(end);
        end
    else
        for v = 1:length(Match2)
            Rate{Match2(v)} = NaN;
        end
    end
end
 %% Convert the Intensity Values to Log Scale
Intensity = cell(1,length(SignalFiles));
 for i=1:length(IdxPostNorm)
    Intensity{IdxPostNorm(i)} = (intensity{IdxPostNorm(i)});
end
%% Build Response Mat 
 AllIntensity = unique(IntVec); RateMat = nan(length(Intensity),length(AllIntensity));
% for r=1:length(IdxPostNorm)
% Ind = find(ismember(AllIntensity,Intensity{IdxPostNorm(r)}));
% RateMat(IdxPostNorm(r),Ind) = Rate{IdxPostNorm(r)};
% end
% Add 0,0 to all units
[0;unique(round(IntVec,2))]; %AllIntensity = [0;unique(IntVec)]; 
RateMat = nan(length(Intensity),length(AllIntensity));
for r=1:length(IdxPostNorm)
Ind = find(ismember(AllIntensity,Intensity{IdxPostNorm(r)}));
RateMat(IdxPostNorm(r),[1;Ind]) = [0;Rate{IdxPostNorm(r)}];
end
%% Plot Dots
close all
for i=1:length(IdxPostNorm)
    fig{i} = figure();
    if answer == '1'
        semilogx(Intensity{IdxPostNorm(i)}*1000,Rate{IdxPostNorm(i)},'*')
        xlim([0, 3250]);
        ylim([0 max(Rate{IdxPostNorm(i)})+0.1]);
        axis square

    else
        plot(LogIntensity{i},Rate{i},'*',Intensity{i},ones(1,length(Intensity{i}))*Spon(i),'--')
    end
end
%% Fit Single Units
% RateAvg{1} = mean([Rate{1};Rate{2}],1); % 4ms
% RateAvg{2} = mean([Rate{3};Rate{4}],1); % 10ms
% RateAvg{3} = mean([Rate{1};Rate{2};Rate{3}(:,1:end-2);Rate{4}(:,1:end-2)]); % all
%IntensityAvg{1} = Intensity{1}; IntensityAvg{2} = Intensity{3}; IntensityAvg{3} = Intensity{1};
% Fit as it should be..
% for i=1:length(fig)
%     figure(fig{i});
%     % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}
%     fun=@(x)sum((Rate{IdxPostNorm(i)}-(x(1)./(1+x(2)*exp(-x(3)*Intensity{IdxPostNorm(i)})))).^2);
%     [x Eval(i)]=fminsearch(fun, [10 10 10]);
%     %options = optimset('MaxFunEvals',10000,'MaxIter',10000);
%     %[x Eval(i)]=fminsearch(fun, [max(Rate{i}) 1 1],options);
%     int_fit{i}=linspace(0,3.250,1000);
%     rate_fit{i} =(x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
%     %rate_fit_norm{i} = (rate_fit{i})/max(rate_fit{i}); % Scale the fitted rate to the maximal response.
%     %SponNorm(i) = Spon(i)/max(rate_fit{i}); % Scale the baseline activity to the maximal response.
%     rate_fit_norm{i} = (rate_fit{i});
%     %SponNorm(i) = Spon(i);
% end

% Forced fits to 0 and 1
for i=1:length(fig)
    figure(fig{i});
    % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}
    fun=@(x)sum((Rate{IdxPostNorm(i)}-(1./(1+1*exp(-x(1)*Intensity{IdxPostNorm(i)})))).^2);
    [x Eval(i)]=fminsearch(fun, 1);
    %options = optimset('MaxFunEvals',10000,'MaxIter',10000);
    %[x Eval(i)]=fminsearch(fun, [max(Rate{i}) 1 1],options);
    int_fit{i}=linspace(0,3.250,1000);
    rate_fit{i} =(1./(1+1*exp(-x(1)*int_fit{i})));
    %rate_fit_norm{i} = (rate_fit{i})/max(rate_fit{i}); % Scale the fitted rate to the maximal response.
    %SponNorm(i) = Spon(i)/max(rate_fit{i}); % Scale the baseline activity to the maximal response.
    rate_fit_norm{i} = (rate_fit{i});
    %SponNorm(i) = Spon(i);
    hold on
    plot(int_fit{i}*1000,rate_fit{i})
    set(gca,'xscale','log')
    xlim([0 max(int_fit{i}*1000)])
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
RateVec = zeros(1,length(Rate)); Limit = 3;
for i = 1:length(Rate)
    if max(Rate{i})<Limit % Select Only Units with small variations in the response
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

AllIntensity2 = round(AllIntensity,2);
[A,B,C] = unique(AllIntensity2);
CountMat = NaN(size(RateMat,1),length(A));
for i=1:length(A)
ind = find(AllIntensity2 == A(i));
for k=1:size(CountMat,1)
CountMat(k,i) = mean(RateMat(k,ind),'omitnan');
end
end
% Clac for Unrounded Int Vals
% for i = 1:length(Durs)  
% AvgCount(i,:) = mean(RateMat(find(StimDur == Durs(i)),:),'omitnan');
% SEMCount(i,:) = std(RateMat(find(StimDur == Durs(i)),:),'omitnan')/sqrt(length(RateMat(find(StimDur == Durs(i)))));
% AvgInt{i} = AllIntensity(find( ~isnan(AvgCount(i,:))))';
% end

% for i = 1:length(Durs)  
% AvgCount(i,:) = mean(CountMat(find(StimDur == Durs(i)),:),'omitnan');
% SEMCount(i,:) = std(CountMat(find(StimDur == Durs(i)),:),'omitnan')/sqrt(length(CountMat(find(StimDur == Durs(i)))));
% AvgInt{i} = A(find( ~isnan(AvgCount(i,:))))';
% end

DurInd = str2double(DurVec);
for i = 1:length(Durs)  
AvgCount(i,:) = mean(CountMat(find(StimDur == Durs(i)),:),'omitnan');
SEMCount(i,:) = std(CountMat(find(DurInd == Durs(i)),:),'omitnan')/sqrt(length(CountMat(find(DurInd == Durs(i)))));
AvgInt{i} = A(find( ~isnan(AvgCount(i,:))))';
end
%% Fit Average Results
%Fit as it should be
for k = 1:size(AvgCount,1)
    fun=@(x)sum((rmmissing(AvgCount(k,:))-(x(1)./(1+x(2)*exp(-x(3)*AvgInt{k})))).^2);
    %fun=@(x)sum((rmmissing(AvgCount(k,:))-(exp(-(AvgInt{k}/x(1)).^x(2))))); % Forced 0
    [x Eval(i)]=fminsearch(fun, [1 1 1]);
    %[x Eval(i)]=fminsearch(fun, [0.001 1000000]); % Forced 0
    AVGint_fit{k}=linspace(min(AvgInt{k}),max(AvgInt{k}),10000);
    AVGrate_fit{k} =(x(1)./(1+x(2)*exp(-x(3)*AVGint_fit{k})));
    %AVGrate_fit{k} =(1-exp(-(AVGint_fit{k}/x(1)).^x(2))); % Forced 0
end

% Forced fits to 0 and 1
% for k = 1:size(AvgCount,1)
%     fun=@(x)sum((rmmissing(AvgCount(k,:))-(1./(1+1*exp(-x(1)*AvgInt{k})))).^2); % Normal sigmoid
%     %fun=@(x)sum((rmmissing(AvgCount(k,:))-(exp(-(AvgInt{k}/0.001).^x(1))))); % Forced 0
%     [x(k) xval(i)]=fminsearch(fun, 1);
%     AVGint_fit{k}=linspace(min(AvgInt{k}),max(AvgInt{k}),10000);
%     AVGrate_fit{k} =(1./(1+1*exp(-x(1)*AVGint_fit{k}))); % Normal sigmoid
%     %AVGrate_fit{k} =(1-exp(-(AVGint_fit{k}/0.001).^x(k))); % Forced 0
% end
%% Plot Average Results
col = ['g','r'];
IdxError = intersect(find(SEMCount(1,:)),find(SEMCount(2,:)));
IdxError2 = IdxError(find(~isnan(SEMCount(1,IdxError))));
[x,y] = intersect(round(AVGint_fit{1},2) , A(IdxError2));
figure();
ax = axes();
for k = 1:length(Durs)
    h{k} = semilogx(AVGint_fit{k}*1000,AVGrate_fit{k},'Color',col(k),'LineWidth',2);
    %legend('AutoUpdate','off')
    hold on
    %[a,b] = intersect(round(AVGint_fit{k},2),round(AllIntensity(IdxError),2),'stable');
    SEMY = AVGrate_fit{k}(y);
    errorbar(x*1000,SEMY,SEMCount(k,IdxError2),'LineStyle','none','Color',col(k))
    %errorbar(AVGint_fit{k}(b)*1000,AVGrate_fit{k}(b),StdCount(k,IdxError)/sqrt(size(RateMat,2)),'.','Color',col(k)) % Plot SEM
    %errorbar(AVGint_fit{k}(b)*1000,AVGrate_fit{k}(b),StdCount(k,IdxError),'.','Color',col(k)) % Plot STD
    %plot(AvgInt{k}*1000,rmmissing(AvgCount(k,:)),'*','Color',col(k));
    hold on
end
legend([h{1},h{2}],{'4ms';'10ms'},"Orientation","vertical","Location","southeast","Box","off");
%xlim([10 max(AVGint_fit{k}*1000)])
%ylim([0 1.1])
xlabel('Intensity [\muW/mm^2]')
ylabel('Scaled Spike Count')
xlim([1 3200]); ylim([0 1.5]);
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
title('Prosthetic Vision')
axes(ax)
% Plot Dots on Average Fits
% for i=1:length(answer2) % loop over all durs
%     a = str2num(answer2{i}); % Convert relevant unit idxs in current stim dur to nums
%     for k = 1:length(a)
%         hold on
%         if i==1
%             semilogx([0;Intensity{a(k)}]*1000,rmmissing(RateMat(a(k),:)),'*g')
%         else
%             semilogx([0;Intensity{a(k)}]*1000,rmmissing(RateMat(a(k),:)),'*r')
%         end
%     end
% end
% legend([h(1) h(2)],['4ms ';'10ms'],'Location','southeast',"box","off")
% Plot Std on Average Fits
%%
save('C:\Users\Itay\Desktop\Yossi Mandel Lab\Extracellular Recordings\Results\Data\Prosthetic\Post-Fit\Data_10msMaxIntensityScale3.mat')