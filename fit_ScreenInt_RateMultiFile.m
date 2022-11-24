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
    Rate{i} = (rate{end}-rate{1}(1))/max(rate{end}); % Subtract baseline activity from all responses and normalize to maximal response    
    Spon(i) = rate{1}(1)/max(rate{end});
    else
    % Use this section for plots without baseline subtraction
    Rate{i} = rate{end}/max(rate{end}); % Scale the fitted rate to the maximal response.
    Spon(i) = rate{1}(1)/max(rate{end}); % Scale the baseline activity to the maximal response.
    end
end
close all
%%
for i=1:length(SignalFiles)
fig{i} = figure();
if answer == '1'
plot(Intensity{i},Rate{i},'*')
else
plot(Intensity{i},Rate{i},'*',Intensity{i},ones(1,length(Intensity{i}))*Spon(i),'--')
end
xlabel('Intensity [cd/m^2]','FontSize',15)
ylabel('Firing Rate [Spikes/s]','FontSize',15)
end
%% fit the data to a sigmoid  {\displaystyle f(x)={\frac {1}{1+e^{-x}}}}

for i=1:length(fig)
    figure(fig{i});
    % fit the data to a sigmoid{\displaystyle f(x)={\frac {1}{1+e^{-x}}}}

    fun=@(x)sum((Rate{i}-(x(1)./(1+x(2)*exp(-x(3)*Intensity{i})))).^2);
    [x Eval(i)]=fminsearch(fun, [100 1 1]);
    %int_fit{i}=linspace(Intensity{i}(1),42.7,100); % Fit with minimal value presented
    int_fit{i}=linspace(0,42.7,100); % Fit with 0 value sa minimum
    rate_fit{i} =(x(1)./(1+x(2)*exp(-x(3)*int_fit{i})));
    rate_fit_norm{i} = (rate_fit{i});
    if answer =='0'
    SponNorm(i) = Spon(i);
    end
    hold on
    plot(int_fit{i},rate_fit{i})

end
%% Find Stimulation Threshold
Thresh = StimThresh(SponNorm);
%% Plot Scaled Fits of All Units
for i=1:length(fig)
count=1;
    for k=1:9:100
    ScaledRate(i,count) = rate_fit_norm{i}(1,k);
    Int(count) = int_fit{i}(1,k);
    count=count+1;
    end
end
ScaledStd = std(ScaledRate)/sqrt(length(ScaledRate));
ScaledAvg = mean(ScaledRate);
figure();
errorbar(Int,ScaledAvg,ScaledStd)
hold on
ScaledSpon = ones(1,length(ScaledAvg))*Thresh;
plot(Int,ScaledSpon,'r--')
ylim([0 1.1]);
xlabel('Intensity [cd/m^2]','FontSize',20)
ylabel('Firing Rate[Spikes/s]','FontSize',20)