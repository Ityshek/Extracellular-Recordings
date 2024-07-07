clear all
close all
% Define path for data
Dir=uigetdir('*.fig','Select a Folder to Load Data Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = [Dir,'\'];

% Define variables
AllInt = []; 

% Iterate over data to build intensity Vector
for i= 1:length(SignalFiles)
    load([pathname,SignalFiles(i).name])
    AllInt = [AllInt;Data.ProstheticIntensity];
    %AllInt = [AllInt;Data.NaturalIntensity];
    close all
end
UAllInt = sort(unique(AllInt));

% Iterate over data to build threshold matrix
count4 = 0; count10 = 0; count2 = 0;DurVec = zeros(1,12);
for i=1:length(SignalFiles)
if ~isempty(regexp(SignalFiles(i).name,'4ms')) == 1 
    count4 = count4+1; DurVec(i) = 4;
elseif ~isempty(regexp(SignalFiles(i).name,'10ms')) == 1
    count10 = count10+1; DurVec(i) = 10;
else
count2 = count2+1; DurVec(i) = 2;
end
end
ThreshMat4 = nan(length(SignalFiles),length(UAllInt)); ThreshMat10 = nan(length(SignalFiles),length(UAllInt));
ThreshMat2 = nan(length(SignalFiles),length(UAllInt));
for i=1:length(SignalFiles)
    load([pathname,SignalFiles(i).name])
    close all
    ActiveIntsIdx = min(find(Data.StimThresh{1}));
    Idx = find(Data.ProstheticIntensity(ActiveIntsIdx) == UAllInt);
    %Idx = find(Data.NaturalIntensity(ActiveIntsIdx) == UAllInt);
    if DurVec(i) == 10
        ThreshMat10(i,Idx:end) = 1; ThreshMat10(i,1:Idx-1) = 0;
    elseif DurVec(i) == 4
        ThreshMat4(i,Idx:end) = 1; ThreshMat4(i,1:Idx-1) = 0;
    else
        ThreshMat2(i,Idx:end) = 1; ThreshMat2(i,1:Idx-1) = 0;
    end
end
ThreshVec{1} = mean(ThreshMat10,1,"omitnan"); ThreshVec{2} = mean(ThreshMat4,1,"omitnan"); ThreshVec{3} = mean(ThreshMat2,1,"omitnan");
%ThreshMat = nan(2,length(SignalFiles),length(UAllInt));
% for i= 1:length(SignalFiles)
%     load([pathname,SignalFiles(i).name])
%     close all
%     if Data.StimDur == '4'
%     IntIdx = find(ismember(UAllInt,Data.ProstheticIntensity));
%     ThreshMat(1,i,IntIdx) = Data.StimThresh{1,1};
%     elseif Data.StimDur == '10'
%     IntIdx = find(ismember(UAllInt,Data.ProstheticIntensity));
%     ThreshMat(2,i,IntIdx) = Data.StimThresh{1,1};
%     else
%     
%     
%     end
%      StrDur(i) = Data.ProstheticIntensity(min(find(Data.StimThresh{1,1})));
%      DurIndx(i) = str2num(Data.StimDur); 
% end
%A(:,:) = mean(ThreshMat,2,'omitnan');


%IntforFit{1} = UAllInt(find(~isnan(A(1,:)))); IntforFit{2} = UAllInt(find(~isnan(A(2,:))));
%% Fit the average thresh vectors
f = figure;
ax = axes();
col = ['r','r']; dotcolor = ['k','c']; % Prosthetic
%col = ['g','g','g']; dotcolor = ['k','c','y']; % Natural
styles = {'-','--',':'};
for i=1:size(ThreshVec,2)
%     % Fit Values
%     fun=@(x)sum((ThreshVec{i}-(x(1)./(1+x(2)*exp(-x(3)*UAllInt')))).^2);
%     %fun=@(x)sum((rmmissing(A(i,:))-(x(1)./(1+x(2)*exp(-x(3)*IntforFit{i}')))).^2);
%     [x Eval(i)]=fminsearch(fun, [0.1 1 1]);
%     int_fit{i}=linspace(0,max(UAllInt),10000);
%     ThreshFit{i} = (x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));

% Forced fits to 0 and 1
for k = 1:size(ThreshVec,2)
    clear x
    fun=@(x)sum((ThreshVec{i}-(1./(1+x(2)*exp(-x(1)*UAllInt')))).^2); % Normal sigmoid
    %fun=@(x)sum((ThreshVec{i}-(exp(-(UAllInt'/x(1)).^x(2))))); % Forced 0
    [x xval(i)]=fminsearch(fun, [0.1 1.2]);
    int_fit{k}=linspace(0,max(UAllInt'),1000000000);
    ThreshFit{k} =(1./(1+x(2)*exp(-x(1)*int_fit{k}))); % Normal sigmoid
    %ThreshFit{k} =(1-exp(-(AVGint_fit{k}/x(1)).^x(2))); % Forced 0
end


%     ZeroIdx = UAllInt(min(find(ThreshVec{i})));
%     ThreshFit{i}(1:max(find(int_fit{i}<ZeroIdx))) = 0;
    % Look for _% to determine Threshold 
    Val = 0.5;
%    ThreshX(i) = int_fit{i}(min(find(ThreshFit{i}>Val)))*1000;
 %   ThreshY(i) = ThreshFit{i}(min(find(ThreshFit{i}>Val)));
    % Plot
    h{i} = semilogx(int_fit{i}*1000,ThreshFit{i}*100,'Color',col(i),LineWidth=2,LineStyle=styles{i});
    %h{i} = semilogx(int_fit{i},ThreshFit{i}*100,'Color',col(i),LineWidth=2,LineStyle=styles{i});
    hold on
    j{i} = plot(UAllInt*1000,ThreshVec{i}*100,'*','Color',dotcolor(i));
    %j{i} = plot(UAllInt,ThreshVec{i}*100,'*','Color',dotcolor(i));
    hold on
end
ylim([0 101])
legend([h{1} h{2} j{1} j{2}],['10ms';' 4ms';'10ms';' 4ms'],"Orientation","vertical","Location","southeast","Box","off")
%legend([h{1} h{2} h{3} j{1} j{2} j{3}],['10ms';' 4ms';' 2ms';'10ms';' 4ms';' 2ms'],"Orientation","vertical","Location","southeast","Box","off")
xlabel('Intensity [\muW/mm^2]','FontSize',20)
ylabel('Percentage of Activated Cells','FontSize',20)
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
%ax.XLim = [1,1000]; 
ax.XLim = [1,3500]; 
axes(ax)
%% Dynamic Range Calc
for i=1:length(ThreshFit)
IndMin(i) = min(find(ThreshFit{i} >= 0.1));
IndMax(i) = min(find(ThreshFit{i} >= 0.9));

ValDMin(i) =  int_fit{i}(IndMin(i))*1000;
ValDMax(i) =  int_fit{i}(IndMax(i))*1000;

DynamicRange(i) =  20*log10(ValDMax(i))/ValDMin(i); %Calculate Dynamic Range in dB.
% bitDepth(i) = log(DynamicRange(i));
end
%% Strength Duration Curve
Durs = unique(DurIndx);
figure(); ax = axes();
for i = 1:length(Durs)
    ThreshVals{i} =  StrDur(find(DurIndx == Durs(i)));
    ThreshValsAvg(i) = mean(ThreshVals{i}*1000);
    ThreshValsStd(i) =std(ThreshVals{i}*1000);
    scatter(DurIndx,StrDur*1000)
    hold on
end
bar(Durs,ThreshValsAvg,'FaceColor','none')
errorbar(Durs,ThreshValsAvg,ThreshValsStd)
xlabel('Duration [ms]') 
ylabel('Intensity [\muW/mm^2]')
xticks(Durs)
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
axes(ax)