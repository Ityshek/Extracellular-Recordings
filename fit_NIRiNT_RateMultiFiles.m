clear all
close all
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.fig'));
pathname = Dir;
Scale = [0,0.00042,0.02,0.27,0.3,0.9,1.01,1.45,1.47,1.93,2.09,2.38,2.8,3.03,3.15,3.6,4.05];
RateMat{1} = zeros(length(SignalFiles),length(Scale)); RateMat{2} = zeros(length(SignalFiles),length(Scale));
RateMatScale{1} = zeros(length(SignalFiles),length(Scale)); RateMatScale{2} = zeros(length(SignalFiles),length(Scale));
UnitCount = zeros(1,2);
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
    Intensity{i} = intensity{1};
    rate=get(h,'ydata');
    rawrate{i} = rate{2};
    Rate{i} = rate{2}/max(rate{2});
    Spon(i) = rate{1}(1)/max(rate{2});
    %Rate{i} = [rate{1}(1),Rate{i}]; Intensity{i} = [0,Intensity{i}]; % Add Baseline Activity as response at 0.
    
% Substracts Baseline Activity of Every Unit, Zero Any Negative Values. 
    % Rate{i} = Rate{i} - rate{1}(1);    
%     for z = 1:length(Rate{i})
%         if Rate{i}(z) < 0
%             Rate{i}(z) = 0;
%         end
%     end
    
    UnitIdx(i) = Unit;
    pause();
end
%%
close all
for i=1:length(SignalFiles)
fig{i} = figure();
plot(Intensity{i},Rate{i},'*',Intensity{i},ones(1,length(Intensity{i}))*Spon(i),'--')
end
%%
for i=1:length(fig)
    figure(fig{i});
    % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}

    fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*Intensity{i})))).^2);
    [x Eval(i)]=fminsearch(fun, [100 1 2]);
    int_fit{i}=linspace(0,4.05,500);
    rate_fit{i} =(x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
    %rate_fit_norm{i} = (rate_fit{i})/max(rate_fit{i}); % Scale the fitted rate to the maximal response.
    %SponNorm(i) = Spon(i)/max(rate_fit{i}); % Scale the baseline activity to the maximal response. 
    rate_fit_norm{i} = (rate_fit{i});
    SponNorm(i) = Spon(i);
    hold on
    plot(int_fit{i},rate_fit{i})
end
%% Find Stimulation Threshold
for k=1:2
thresh(k) = StimThresh(SponNorm(UnitIdx==k));
UnitTresh = cell(1,2);
Thresh = mean(thresh);
for i=1:length(fig)
UnitTresh{UnitIdx(i)} = [UnitTresh{UnitIdx(i)};int_fit{i}(min(find(rate_fit_norm{i}>thresh(k))))];
end
end
PulseDurSEM(1) = std(UnitTresh{1})/sqrt(length(UnitTresh{1}));
PulseDurSEM(2) = std(UnitTresh{2})/sqrt(length(UnitTresh{2}));
%%
for i =  1:length(fig)
Actthresh(i) = int_fit{i}(min(find(rate_fit_norm{i}>=SponNorm(i))));
end
%%
color = ['r','b'];
UnitType = {'-','--'};
prompt = {'Select Units to include in further analysis'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'[1:1:15]'};
answer = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
Thresh = mean(SponNorm)

RatefitUnitMat = cell(1,2); RatefitNormUnitMat = cell(1,2);
figure();
for i=1:length(answer) 
subplot(1,2,1)
RatefitMat(i,:) = rate_fit{answer(i)};
RatefitUnitMat{UnitIdx(answer(i))} =[RatefitUnitMat{UnitIdx(answer(i))};rate_fit{answer(i)}]; 
plot(int_fit{i},rate_fit{answer(i)},[UnitType{UnitIdx(answer(i))},color(UnitIdx(i))]...
    ,'DisplayName',['Fit',num2str(answer(i)),' Eval=',num2str(Eval(answer(i)))])
xlabel('Intensity [mW/mm^2]','FontSize',15)
ylabel('Firing Rate[Hz]','FontSize',15)
title('Single Units Raw Fitted');
hold on
subplot(1,2,2)
RatefitNormMat(i,:) = rate_fit_norm{answer(i)};
RatefitNormUnitMat{UnitIdx(answer(i))} =[RatefitNormUnitMat{UnitIdx(answer(i))};rate_fit_norm{answer(i)}]; 
plot(int_fit{i},rate_fit_norm{answer(i)},[UnitType{UnitIdx(answer(i))},color(UnitIdx(i))]...
    ,'DisplayName',['Fit',num2str(answer(i)),' Eval=',num2str(Eval(answer(i)))])
ylim([0 1.1])
title('Single Units Scaled Fitted');
hold on
end
legend();
xlabel('Intensity[mW/mm^2]','FontSize',15)
ylabel('Firing Rate[Hz]','FontSize',15)

% figure();
% %plot(int_fit{1},mean(RatefitMat,1),'LineWidth',2) % Plot Total average response (fitted)
% hold on
% plot(int_fit{1},ones(1,length(int_fit{1}))*Thresh,'r--') % Plot Total average  baseline activity
% hold on
% plot(int_fit{1},mean(RatefitUnitMat{1},1),'b-',int_fit{1},mean(RatefitUnitMat{2},1),'b--')
% legend('Baseline Activity','10ms','4ms')
% xlabel('Intensity[mW/mm^2]','FontSize',15)
% ylabel('Firing Rate[Hz]')
figure();
%plot(int_fit{1},mean(RatefitNormMat,1),'LineWidth',2) % Plot Total average response (fitted & scaled)
hold on
plot(int_fit{1},ones(1,length(int_fit{1}))*Thresh,'r--') % Plot Total average scaled baseline activity
ylim([0.1 1.1])
hold on
plot(int_fit{1},mean(RatefitNormUnitMat{1},1),'b-',int_fit{1},mean(RatefitNormUnitMat{2},1),'b--')
xlabel('Intensity[mW/mm^2]','FontSize',15)
ylabel('Firing Rate[Hz]')
PulseDurThresh(1) = int_fit{i}(min(find(mean(RatefitNormUnitMat{1},1)>thresh(1))));
PulseDurThresh(2) = int_fit{i}(min(find(mean(RatefitNormUnitMat{2},1)>thresh(2))));
%%
for k=1:length(RatefitNormUnitMat)
SEMNorm = std(RatefitNormUnitMat{k})/sqrt(size(RatefitNormUnitMat{k},2));
points = round(linspace(0.01,3,7),1);
for p=1:length(points)
    pointidx{k}(p) = min(find(round(int_fit{1},1) == points(p)));
    SEMPoints{k}(p) = SEMNorm(pointidx{k}(p));
end
end
%% Draw SEM values at different points on the curves
errorbar(int_fit{1}(pointidx{1}),mean(RatefitNormUnitMat{1}(:,pointidx{1})),SEMPoints{1},'b.')
hold on
errorbar(int_fit{1}(pointidx{2}),mean(RatefitNormUnitMat{2}(:,pointidx{2})),SEMPoints{2},'b.')
hold on
plot(int_fit{1},ones(1,length(int_fit{1}))*Thresh,'r--')
ylim([0 1.1])
xlim([0 3])
legend('Stimulation Threshold','10ms','4ms')