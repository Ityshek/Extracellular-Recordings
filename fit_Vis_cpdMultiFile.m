close all
clear all
Dir=uigetdir('*.fig','Select a Folder to Load Fig Files From');
SignalFiles=dir(fullfile(Dir, '*.fig'));
pathname = Dir;

for i=1:length(SignalFiles)
    openfig([pathname,'\',SignalFiles(i).name ]);
    h=get(gca,'children');
    cpd=get(h,'xdata');
    CPD{i} = cpd{1};
    rate=get(h,'ydata');
    Rate{i} = rate{end};
    Spon(i) = rate{1}(1);
pause();
end

%%
close all
for i=1:length(SignalFiles)
fig{i} = figure();
plot(CPD{i},Rate{i},'*')
xlabel('CPD','FontSize',20)
ylabel('Firing Rate [Spikes/s]','FontSize',20)
ylim([0 70])
end
%% fit the data to a sigmoid  {\displaystyle \Gamma _{\sigma _{1},\sigma _{2}}(x)=I*{\frac {1}{\sigma _{1}{\sqrt {2\pi }}}}\,e^{-{\frac {x^{2}}{2\sigma _{1}^{2}}}}-I*{\frac {1}{\sigma _{2}{\sqrt {2\pi }}}}\,e^{-{\frac {x^{2}}{2\sigma _{2}^{2}}}}.}
% fun=@(x)sum((rate-(x(1)*(1/((2*pi)^0.5*x(2)))*exp((-CPD.^2-x(4))./(2*(x(2)).^2))-x(6)*(1/((2*pi)^0.5*x(3)))*exp((-CPD.^2-x(5))./(2*(x(3)).^2)))).^2);
for i=1:length(fig)
    figure(fig{i});
fun=@(x)sum((Rate{i}-(x(1)*(1/((2*pi)^0.5*x(2)))*exp(-(CPD{i}-x(3)).^2./(2*(x(2)).^2))-x(4)*(1/((2*pi)^0.5*x(5)))*exp(-(CPD{i}-x(6)).^2./(2*(x(5)).^2)))).^2);

 %[x eval(i)]=fminsearch(fun, [1 0.2 1 10 1 1]);
 [x Eval(i)]=fminsearch(fun, [5 0.5 1 -1 1 3]);
cpd_fit{i}=linspace(CPD{i}(1),0.6,100);
% rate_fit =(x(1)*(1/((2*pi)^0.5*x(2)))*exp((-cpd_fit.^2-x(4))./(2*(x(2)).^2))-x(6)*(1/((2*pi)^0.5*x(3)))*exp((-cpd_fit.^2-x(5))./(2*(x(3)).^2)));
rate_fit{i}=(x(1)*(1/((2*pi)^0.5*x(2)))*exp(-(cpd_fit{i}-x(3)).^2./(2*(x(2)).^2))-x(4)*(1/((2*pi)^0.5*x(5)))*exp(-(cpd_fit{i}-x(6)).^2./(2*(x(5)).^2)));
rate_fit_norm{i} = (rate_fit{i})/max(rate_fit{i}); % Scale the fitted rate to the maximal response.
SponNorm(i) = Spon(i)/max(rate_fit{i}); % Scale the baseline activity to the maximal response.
hold on
plot(cpd_fit{i},rate_fit{i},cpd_fit{i},ones(1,length(cpd_fit{i}))*Spon(i),'--')
end
%% Plot Scaled Fits of All Units
figure();
for i=1:length(fig)
plot(cpd_fit{1},rate_fit_norm{i})
hold on
end
plot(cpd_fit{1},ones(1,length(cpd_fit{i}))*mean(SponNorm),'r--')
ylim([0 1.1])
xlabel('CPD','FontSize',20)
ylabel('Scaled Firing Rate','FontSize',20)
