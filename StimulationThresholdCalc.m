clear all
close all
% Define path for data
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = [Dir,'\'];

% Define variables
AllInt = []; 

% Iterate over data to build intensity Vector
for i= 1:length(SignalFiles)
    load([pathname,SignalFiles(i).name])
    AllInt = [AllInt;Data.ProstheticIntensity];
    close all
end
UAllInt = sort(unique(AllInt));

% Iterate over data to build threshold matrix
ThreshMat = nan(2,length(SignalFiles),length(UAllInt));
for i= 1:length(SignalFiles)
    load([pathname,SignalFiles(i).name])
    close all
    if Data.StimDur == '4'
    IntIdx = find(ismember(UAllInt,Data.ProstheticIntensity));
    ThreshMat(1,i,IntIdx) = Data.StimThresh{1,1};
    elseif Data.StimDur == '10'
    IntIdx = find(ismember(UAllInt,Data.ProstheticIntensity));
    ThreshMat(2,i,IntIdx) = Data.StimThresh{1,1};
    end
end
A(:,:) = mean(ThreshMat,2,'omitnan');

IntforFit{1} = UAllInt(find(~isnan(A(1,:)))); IntforFit{2} = UAllInt(find(~isnan(A(2,:))));
%% Fit the average thresh vectors
f = figure;
ax = axes();
col = ['r','b'];
for i=1:size(A,1)
    % Fit Values
    fun=@(x)sum((rmmissing(A(i,:))-(x(1)./(1+x(2)*exp(-x(3)*IntforFit{i}')))).^2);
    [x Eval(i)]=fminsearch(fun, [1 1 1]);
    int_fit{i}=linspace(0,max(UAllInt),5000);
    ThreshFit{i} = (x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
    % Look for _% to determine Threshold 
    Val = 0.5;
    ThreshX(i) = int_fit{i}(min(find(ThreshFit{i}>Val)))*1000;
    ThreshY(i) = ThreshFit{i}(min(find(ThreshFit{i}>Val)));
    % Plot
    h{i} = semilogx(int_fit{i}*1000,ThreshFit{i},'Color',col(i));
    hold on
    plot(ThreshX(i),ThreshY(i),'*','Color','k')
    hold on
end
ylim([0 1.1])
legend([h{1} h{2}],['4ms ';'10ms'],Location='southeast',Orientation='horizontal')
xlabel('Intensity [\muW/mm^2]','FontSize',20)
ylabel('Activation Probability','FontSize',20)
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.XLim = [10,3500];
axes(ax)