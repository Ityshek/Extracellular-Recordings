clear all
clc
close all
[filename, pathname]=uigetfile();
openfig([pathname filename ]);

h=get(gca,'children');
CPD=get(h,'xdata');
rate=get(h,'ydata');
figure; 
plot(CPD,rate,'*r')

%% fit the data to a sigmoid  {\displaystyle \Gamma _{\sigma _{1},\sigma _{2}}(x)=I*{\frac {1}{\sigma _{1}{\sqrt {2\pi }}}}\,e^{-{\frac {x^{2}}{2\sigma _{1}^{2}}}}-I*{\frac {1}{\sigma _{2}{\sqrt {2\pi }}}}\,e^{-{\frac {x^{2}}{2\sigma _{2}^{2}}}}.}
% fun=@(x)sum((rate-(x(1)*(1/((2*pi)^0.5*x(2)))*exp((-CPD.^2-x(4))./(2*(x(2)).^2))-x(6)*(1/((2*pi)^0.5*x(3)))*exp((-CPD.^2-x(5))./(2*(x(3)).^2)))).^2);
fun=@(x)sum((rate-(x(1)*(1/((2*pi)^0.5*x(2)))*exp(-(CPD-x(3)).^2./(2*(x(2)).^2))-x(4)*(1/((2*pi)^0.5*x(5)))*exp(-(CPD-x(6)).^2./(2*(x(5)).^2))   )).^2);

% [x eval]=fminsearch(fun, [1 0.2 1 10 1 1]);
[x eval]=fminsearch(fun, [5 0.2 1 -1 1 3  ]);

cpd_fit=linspace(0.1*CPD(1),CPD(end),100);
% rate_fit =(x(1)*(1/((2*pi)^0.5*x(2)))*exp((-cpd_fit.^2-x(4))./(2*(x(2)).^2))-x(6)*(1/((2*pi)^0.5*x(3)))*exp((-cpd_fit.^2-x(5))./(2*(x(3)).^2)));
rate_fit=(x(1)*(1/((2*pi)^0.5*x(2)))*exp(-(cpd_fit-x(3)).^2./(2*(x(2)).^2))-x(4)*(1/((2*pi)^0.5*x(5)))*exp(-(cpd_fit-x(6)).^2./(2*(x(5)).^2)));
hold on 
plot(cpd_fit,rate_fit)