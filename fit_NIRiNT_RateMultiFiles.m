clear all
close all

Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.fig'));
pathname = Dir;
Scale = [0.00042,0.02,0.27,0.3,0.9,1.01,1.45,1.47,1.93,2.09,2.38,2.8,3.03,3.15,3.6,4.05]*1000; % Intensity Vals in micro-Watts/mm^2
RateMat{1} = zeros(length(SignalFiles),length(Scale)); RateMat{2} = zeros(length(SignalFiles),length(Scale));
RateMatScale{1} = zeros(length(SignalFiles),length(Scale)); RateMatScale{2} = zeros(length(SignalFiles),length(Scale));
UnitCount = zeros(1,2);

prompt = {'1 = Baseline Activity subtraction | 0 = Show Baseline Activity'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1'};
answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));

MinInt = 0.5; MaxInt = 1;
AllInt = [];     ResponseMat = NaN(length(SignalFiles),21);

for i=1:length(SignalFiles)
    rate = {};
    if ~isempty(regexp(SignalFiles(i).name,'10ms'))
        Unit = 1;
        RateMat{2}(i,:) = nan; RateMatScale{2}(i,:) = nan;
        UnitCount(1) = UnitCount(1) + 1;
    else
        Unit = 2;
        RateMat{1}(i,:) = nan; RateMatScale{1}(i,:) = nan;
        UnitCount(2) = UnitCount(2) + 1;
    end



    openfig([pathname,'\',SignalFiles(i).name ]);
    %fig{i} = ans;
    %title([SignalFiles(i).name]);
    h=get(gca,'children');
    intensity=get(h,'xdata');    
    %Intensity{i} = intensity{1}*1000; %Convert Vals to micro-Watts/mm^2
    
    Intensity{i} = intensity{1};    
    if min(Intensity{i})<MinInt
        MinInt = min(Intensity{i});
    end
    if max(Intensity{i}) > MaxInt
        MaxInt = max(Intensity{i});
    end
    AllInt = [AllInt,Intensity{i}];
    UAllInt = unique(AllInt);
    


    rate=get(h,'ydata');
       if answer == '1'
    ResponseRate{i} = (rate{2}-rate{1}(1)); %Substracts Baseline Activity of Every Unit
       end


    %     if Intensity{i}(1) == 0
    %         Intensity{i} = Intensity{i}(2:end);
    %         ResponseRate{i} = ResponseRate{i}(2:end);
    %     end
    
    % Zero any negative Vals after baseline subtraction
    for  z = 1:length(ResponseRate{i})


        if ResponseRate{i}(z) <0
            ResponseRate{i}(z) = 0;
        end
    end

 
        Rate{i} = ResponseRate{i}/max(ResponseRate{i}); %  Normalize to Maximal Response.
        %Rate{i} = ResponseRate{i};
%     else
%         Rate{i} = rate{2}/max(rate{2});
%         Spon(i) = rate{1}(1)/max(rate{2});
    



    UnitIdx(i) = Unit;
    pause();
end

for i = 1:length(Intensity)
    for k=1:length(Intensity{i})
        A = find(UAllInt == Intensity{i}(k));
        ResponseMat(i,A) = Rate{i}(k);
    end
end
 %% Convert the Intensity Values to Log Scale
% for i=1:length(Intensity)
%     LogIntensity{i} = log10(Intensity{i});
% end
%%
close all
for i=1:length(SignalFiles)
    fig{i} = figure();
    if answer == '1'
        plot(Intensity{i}*1000,Rate{i},'*')
        %xlim([0, max(Intensity{i})+1]);
    else
        plot(LogIntensity{i},Rate{i},'*',Intensity{i},ones(1,length(Intensity{i}))*Spon(i),'--')
    end
end
%%


% RateAvg{1} = mean([Rate{1};Rate{2}],1); % 4ms
% RateAvg{2} = mean([Rate{3};Rate{4}],1); % 10ms
% RateAvg{3} = mean([Rate{1};Rate{2};Rate{3}(:,1:end-2);Rate{4}(:,1:end-2)]); % all

%IntensityAvg{1} = Intensity{1}; IntensityAvg{2} = Intensity{3}; IntensityAvg{3} = Intensity{1};
for i=1:length(fig)
    figure(fig{i});
    % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}
    fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*Intensity{i})))).^2);
    [x Eval(i)]=fminsearch(fun, [1 1 1]);
    %options = optimset('MaxFunEvals',10000,'MaxIter',10000);
    %[x Eval(i)]=fminsearch(fun, [max(Rate{i}) 1 1],options);
    int_fit{i}=linspace(MinInt,MaxInt,100);
    rate_fit{i} =(x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
    %rate_fit_norm{i} = (rate_fit{i})/max(rate_fit{i}); % Scale the fitted rate to the maximal response.
    %SponNorm(i) = Spon(i)/max(rate_fit{i}); % Scale the baseline activity to the maximal response.
    rate_fit_norm{i} = (rate_fit{i});
    %SponNorm(i) = Spon(i);
    hold on
    plot(int_fit{i}*1000,rate_fit{i})
    set(gca,'xscale','log')
    xlim([1 max(int_fit{i}*1000)])
    xlabel('Intensity [\muW/mm^2]')
    ylabel('Scaled Firing Rate')
end
%%
prompt = {'Select Units to include in further analysis'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'[1,2,3,4,5,6,7]'};
answer1 = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));

for k=1:length(answer1)
DurIdx = strfind(SignalFiles(answer1(k)).name,'ms');
StimDur(k) = str2num(SignalFiles(answer1(k)).name(DurIdx-2:DurIdx-1));
end


Durs = unique(StimDur); % Count num of different durations
CountMat = cell(2,1);
for k=1:length(Durs)
for i = 1:length(answer1)
    if StimDur(i) == Durs(k)
    CountMat{k} = [CountMat{k};rate_fit{answer1(k)}];   
    end
    end
end

% Average across Stim durs in all units and extract relevant Intensities
prompt = {'Insert Idxs for 4ms:','Insert Idxs for 10ms:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1,2,3','5,6'};
answer2 = inputdlg(prompt,dlgtitle,dims,definput);

for i=1:length(answer2) 
AvgCount(i,:) = mean(ResponseMat(min(str2num(answer2{i})):max(str2num(answer2{i})),:),'omitnan');
StdCount(i,:) = std(ResponseMat(min(str2num(answer2{i})):max(str2num(answer2{i})),:),'omitnan');
AvgInt{i} = UAllInt(find( ~isnan(AvgCount(i,:)))); 
end

% Fit Average Results
for k = 1:size(AvgCount,1)
    fun=@(x)sum((rmmissing(AvgCount(k,:))-(x(1)./(1+x(2)*exp(-x(3)*AvgInt{k})))).^2);
    [x Eval(i)]=fminsearch(fun, [1 1 1]);
    AVGint_fit{k}=linspace(MinInt,MaxInt,100);
    AVGrate_fit{k} =(x(1)./(1+x(2)*exp(-x(3)*AVGint_fit{k})));
end

%%
% Plot Average Results
col = ['r','b'];
figure();
ax = axes();
for k = 1:length(Durs)
h{k} = semilogx(AVGint_fit{k}*1000,AVGrate_fit{k},'Color',col(k));
%legend('AutoUpdate','off')
hold on
plot(AvgInt{k}*1000,rmmissing(AvgCount(k,:)),'*','Color',col(k));
end
legend([h{1} h{2}],['4ms ';'10ms'])
xlim([10 max(AVGint_fit{k}*1000)])
ylim([0 1.1])
xlabel('Intensity [\muW/mm^2]')
ylabel('Scaled Spike Count')
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.XLim = [10,3500];
axes(ax)



















%%
% % Fit for AVG
% for i=1:length(RateAvg)
%     fun=@(x)sum((RateAvg{i}-(x(1)./(1+x(2)*exp(-x(3)*IntensityAvg{i})))).^2);
%     options = optimset('MaxFunEvals',10000,'MaxIter',10000);
%     [x Eval(i)]=fminsearch(fun, [max(Rate{i}) 1 1],options);
%     %int_fit{i}=linspace(min(Intensity{answer(i)}),max(Intensity{answer(i)}),500);
%     rate_fit_avg{i} =(x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
% end
% figure();
% plot(int_fit{1},rate_fit_avg{1},'Color','r')
% hold on
% plot(int_fit{1},rate_fit_avg{2},'Color','b')
% hold on
% plot(int_fit{1},rate_fit_avg{3},'Color','k','LineWidth',2)
% legend('4ms','10ms','Avarage')
% set(gca,'XScale','log')
% xlabel('Intensity Log [\muW/mm^2]','FontSize',15)
% ylabel('Scaled Firing Rate','FontSize',15)
%% Find Stimulation Threshold
% for k=1:2
%     thresh(k) = StimThresh(SponNorm(UnitIdx==k));
%     UnitTresh = cell(1,2);
%     Thresh = mean(thresh);
%     for i=1:length(fig)
%         UnitTresh{UnitIdx(i)} = [UnitTresh{UnitIdx(i)};int_fit{i}(min(find(rate_fit_norm{i}>thresh(k))))];
%     end
% end
% PulseDurSEM(1) = std(UnitTresh{1})/sqrt(length(UnitTresh{1}));
% PulseDurSEM(2) = std(UnitTresh{2})/sqrt(length(UnitTresh{2}));
%%
% for i =  1:length(fig)
%     Actthresh(i) = int_fit{i}(min(find(rate_fit_norm{i}>=SponNorm(i))));
% end
% %%
% 
% color = ['r','b'];
% UnitType = {'-','--'};
% % prompt = {'Select Units to include in further analysis'};
% % dlgtitle = 'Input';
% % dims = [1 35];
% % definput = {'[3,5,13,14]'};
% answer = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
% %Thresh = mean(SponNorm)
% 
% RatefitUnitMat = cell(1,2); RatefitNormUnitMat = cell(1,2);
% figure();
% for i=1:length(answer)
%     subplot(1,2,1)
%     RatefitMat(i,:) = rate_fit{answer(i)};
%     RatefitUnitMat{UnitIdx(answer(i))} =[RatefitUnitMat{UnitIdx(answer(i))};rate_fit{answer(i)}];
%     plot(int_fit{i},rate_fit{answer(i)},[UnitType{UnitIdx(answer(i))},color(UnitIdx(i))]...
%         ,'DisplayName',['Fit',num2str(answer(i)),' Eval=',num2str(Eval(answer(i)))])
%     xlabel('Intensity [mW/mm^2]','FontSize',15)
%     ylabel('Firing Rate[Hz]','FontSize',15)
%     title('Single Units Raw Fitted');
%     hold on
%     subplot(1,2,2)
%     RatefitNormMat(i,:) = rate_fit_norm{answer(i)};
%     RatefitNormUnitMat{UnitIdx(answer(i))} =[RatefitNormUnitMat{UnitIdx(answer(i))};rate_fit_norm{answer(i)}];
%     plot(int_fit{i},rate_fit_norm{answer(i)},[UnitType{UnitIdx(answer(i))},color(UnitIdx(i))]...
%         ,'DisplayName',['Fit',num2str(answer(i)),' Eval=',num2str(Eval(answer(i)))])
%     ylim([0 1.1])
%     title('Single Units Scaled Fitted');
%     hold on
% end
% legend();
% xlabel('Intensity[mW/mm^2]','FontSize',15)
% ylabel('Firing Rate[Hz]','FontSize',15)
% %%
% % figure();
% % %plot(int_fit{1},mean(RatefitMat,1),'LineWidth',2) % Plot Total average response (fitted)
% % hold on
% % plot(int_fit{1},ones(1,length(int_fit{1}))*Thresh,'r--') % Plot Total average  baseline activity
% % hold on
% % plot(int_fit{1},mean(RatefitUnitMat{1},1),'b-',int_fit{1},mean(RatefitUnitMat{2},1),'b--')
% % legend('Baseline Activity','10ms','4ms')
% % xlabel('Intensity[mW/mm^2]','FontSize',15)
% % ylabel('Firing Rate[Hz]')
% %%
% figure();
% semilogx(int_fit{1},mean(RatefitNormMat,1),'LineWidth',2,'Color','r') % Plot Total average response (fitted & scaled)
% set(gca,'xscale','log')
% hold on
% [C, ia, ic] = unique(StimDur);
% Rate4ms = mean(RatefitNormMat(1:2,:),1); Rate10ms = mean(RatefitNormMat(3:4,:),1);
% semilogx(int_fit{1},Rate4ms,'LineWidth',2,'Color','b')
% hold on
% semilogx(int_fit{1},Rate10ms,'LineWidth',2,'Color','k')
% xlabel('Intensity[\muW/mm^2]','FontSize',15)
% ylabel('Firing Rate[Hz]')
% legend('Avarege',num2str(C(1)),num2str(C(2)))
%%
%plot(int_fit{1},ones(1,length(int_fit{1}))*Thresh,'r--') % Plot Total average scaled baseline activity
%ylim([0.1 1.1])
%hold on
%plot(int_fit{1},mean(RatefitNormUnitMat{1},1),'b-',int_fit{1},mean(RatefitNormUnitMat{2},1),'b--')
%semilogx(int_fit{1},mean(RatefitNormUnitMat{1},1),'b-',int_fit{1},mean(RatefitNormUnitMat{2},1),'b--')

%PulseDurThresh(1) = int_fit{i}(min(find(mean(RatefitNormUnitMat{1},1)>thresh(1))));
%PulseDurThresh(2) = int_fit{i}(min(find(mean(RatefitNormUnitMat{2},1)>thresh(2))));
% %%
% for k=1:length(RatefitNormUnitMat)
%     SEMNorm = std(RatefitNormUnitMat{k})/sqrt(size(RatefitNormUnitMat{k},2));
%     points = round(linspace(0.01,3,7),1);
%     for p=1:length(points)
%         pointidx{k}(p) = min(find(round(int_fit{1},1) == points(p)));
%         SEMPoints{k}(p) = SEMNorm(pointidx{k}(p));
%     end
% end
% %% Draw SEM values at different points on the curves
% errorbar(int_fit{1}(pointidx{1}),mean(RatefitNormUnitMat{1}(:,pointidx{1})),SEMPoints{1},'b.')
% hold on
% errorbar(int_fit{1}(pointidx{2}),mean(RatefitNormUnitMat{2}(:,pointidx{2})),SEMPoints{2},'b.')
% hold on
% plot(int_fit{1},ones(1,length(int_fit{1}))*Thresh,'r--')
% ylim([0 1.1])
% xlim([0 3])
% legend('Stimulation Threshold','10ms','4ms')