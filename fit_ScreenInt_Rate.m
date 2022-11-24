close all
clear all
[filename, pathname]=uigetfile();
openfig([pathname filename ]);

h=get(gca,'children');
Intensity=get(h,'xdata');
rate=get(h,'ydata');
a = rate{2};
figure; 
plot(Intensity{1},rate{2},'*r')

%% fit the data to a sigmoid  {\displaystyle f(x)={\frac {1}{1+e^{-x}}}}
fun=@(x)sum((a-(x(1)./(1+x(2)*exp(-x(3)*Intensity{1})))).^2);
[x eval]=fminsearch(fun, [100 2 1]);
int_fit=linspace(Intensity{1}(1),4.05,100);
rate_fit =(x(1)./(1+x(2)*exp(-x(3)*int_fit)));
hold on 
plot(int_fit,rate_fit)