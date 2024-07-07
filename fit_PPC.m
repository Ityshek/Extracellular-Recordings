% Fit CFF Results
close all; clear all;
%% Load Data Files
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.mat'));
pathname = Dir;

Dates = strings; Channels = strings;
Locs = strings; DurVec = strings;
StimCondition = []; Freqvec = [];
for i=1:length(SignalFiles)
    load([pathname,'\',SignalFiles(i).name ])
    Frequencies{i} = Data.StimFreq; Freqvec = [Freqvec;Data.StimFreq];
    PPCs{i} = Data.PPC;
    if regexp(SignalFiles(i).name,['Natural']) == 1 % Check for Natural / Prosthetic
        StimCondition(i) = 1;
    else
        StimCondition(i) = 2;
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
close all
%% Plot PPCs
NVis = count(num2str(StimCondition),'1'); NNir = count(num2str(StimCondition),'2');
titles = {'Natural';'Prosthetic'}; LineStyles = {'-'; '-.'}; Names = {['Natural(N=',num2str(NVis),')'],['Prosthetic(N=',num2str(NNir),')']}; Colors = {'b','r'};
figure();
for i=1:length(PPCs)
%     PPCFigs{i} = figure(); 
%     ax = axes();
    h{i} = plot(Frequencies{i},PPCs{i},'LineStyle',LineStyles{StimCondition(i)},LineWidth=2,Marker='*',Color=Colors{StimCondition(i)})
    ylabel('PPC Value'); xlabel('Frequencies [Hz]'); 
    hold on
%     ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 15;
%     ax.Box = 'off'; ax.Color = "none";
%     if StimCondition(i) == 1
%         title([titles{1}, ' Date:',SortVar(i).Date,' Loc:',SortVar(i).Loc,' Channel:',SortVar(i).Channel])
%     else
%         title([titles{2}, ' Date:',SortVar(i).Date,' Loc:',SortVar(i).Loc,' Channel:',SortVar(i).Channel])
%     end
end
legend([h{1},h{end}],{'Natural';'Prosthetic'});
xlim([0 65]); ylim([0 1]);
legend(EdgeColor='none',Color='none'); 
axis square; box off;
set(gca,'color','none','FontSize',15)
%% Build Average Vectors
AllFreqs = unique(Freqvec); 
PPCMat = NaN(length(PPCs),length(AllFreqs));
% Build matrix of PPCs
for i=1:length(PPCs)
[a,b] = intersect(AllFreqs,Frequencies{i});
PPCMat(i,b) = PPCs{i};
end
for i=1:2
MeanPPCs{i} = mean(PPCMat(find(StimCondition==i),:),'omitnan'); StdPPCs{i} = std(PPCMat(find(StimCondition==i),:),'omitnan'); 
NewFreqs{i} = AllFreqs(find(~isnan(MeanPPCs{i})));
MeanPPCs{i} = rmmissing(MeanPPCs{i}); SEMPPCs{i} = rmmissing(StdPPCs{i})/sqrt(length(rmmissing(StdPPCs{i})));
end
%% Plot Averages
figure();
for i=1:length(Names)
    errorbar(NewFreqs{i},MeanPPCs{i},SEMPPCs{i},LineStyles{StimCondition(i)},LineWidth=2,DisplayName=Names{i},LineStyle=LineStyles{i},Color=Colors{i})
    %plot(NewFreqs{i},MeanPPCs{i},LineStyles{StimCondition(i)},LineWidth=2,DisplayName=Names{i},Marker='*',LineStyle=LineStyles{i})
    hold on
ylabel('PPC Value'); xlabel('Frequencies [Hz]');
legend(EdgeColor='none',Color='none'); 
xlim([0 65]); ylim([0 1]);
axis square; box off;
set(gca,'color','none','FontSize',15)
end
%% Plot Averages with SEMs as Clouds
figure();
NSEMUp = MeanPPCs{1}+SEMPPCs{1}; NSEMDown = MeanPPCs{1}-SEMPPCs{1};

PSEMUp = MeanPPCs{2}+SEMPPCs{2}; PSEMDown = MeanPPCs{2}-SEMPPCs{2};

inBetween =[NSEMUp, fliplr(NSEMDown)];
inBetween2 =[PSEMUp, fliplr(PSEMDown)];
T = [1,2,4,8,16,32,64];
x2 = [T, fliplr(T)];
col=['b','r'];
h1 = plot(T,MeanPPCs{1},'b',LineWidth=2);
hold on
h2 = plot(T,MeanPPCs{2},'r',LineWidth=2);
h3 = fill(x2,inBetween,col(1));
 set(h3,'EdgeColor','none')
 alpha(0.07) 
hold on 
h4 = fill(x2, inBetween2,col(2));
set(h4,'EdgeColor','none')
 alpha(0.07) 
box off
ylabel('PPC Value'); xlabel('Frequencies [Hz]');
legend([h1,h2],Names,EdgeColor='none',Color='none'); 
xlim([0 65]); ylim([0 1]);
axis square; box off;
set(gca,'color','none','FontSize',15)
%% Fit averages
for i=1:2
    fun = @(x)sum((MeanPPCs{i}' - (x(4) + x(1)/x(2).*((NewFreqs{i}-x(3))/x(2)).^(x(1)-1).*exp(-((NewFreqs{i}-x(3))./x(2)).^x(1)))).^2);
    options = optimset('MaxFunEvals',10000,'MaxIter',10000);
    [x Eval(i)]=fminsearch(fun, [0.2 4 0.5 0.1],options);
    %[x Eval(i)]=fminsearch(fun, [1 10 4]);
    Freq_fit{i}=linspace(1,64,1000);
    PPC_fit{i} =x(4) + x(1)/x(2).*((Freq_fit{i}-x(3))/x(2)).^(x(1)-1).*exp(-((Freq_fit{i}-x(3))./x(2)).^x(1));
    
    %[param{i},ci{i}] = wblfit(MeanPPCs{i});
    %pd{i} = fitdist(MeanPPCs{i}','wbl');
end
%% Plot Averages 
figure(); colors = ['r';'b'];
for i=1:2
    plot(NewFreqs{i},MeanPPCs{i},'*',"Color",colors(i))
    hold on
%     plot(pd{i},"PlotType","pdf","Discrete")
    plot(Freq_fit{i},PPC_fit{i},'-','Color',colors(i))
    
end
