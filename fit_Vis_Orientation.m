clear all
clc
close all
[filename, pathname]=uigetfile();
openfig([pathname filename ]);
h=get(gca,'children');
orientation=get(h,'xdata');
rate=get(h,'ydata');
figure; 
plot(orientation,rate,'*r')

%% fit the data to a sigmoid  {\displaystyle \Gamma _{\sigma _{1},\sigma _{2}}(x)=I*{\frac {1}{\sigma _{1}{\sqrt {2\pi }}}}\,e^{-{\frac {x^{2}}{2\sigma _{1}^{2}}}}-I*{\frac {1}{\sigma _{2}{\sqrt {2\pi }}}}\,e^{-{\frac {x^{2}}{2\sigma _{2}^{2}}}}.}
  fun=@(x)sum((rate-(x(1)+x(2)*exp(-((orientation-x(3)).^2)/(2*x(4)^2))+x(5)*exp(-((orientation-x(6)).^2)/(2*x(4)^2)))).^2);

 [x eval]=fminsearch(fun, [10 1 60 10 1 100]);

orien_fit=linspace(orientation(1),orientation(end),100);
% rate_fit =(x(1)*(1/((2*pi)^0.5*x(2)))*exp((-cpd_fit.1 ^2-x(4))./(2*(x(2)).^2))-x(6)*(1/((2*pi)^0.5*x(3)))*exp((-cpd_fit.^2-x(5))./(2*(x(3)).^2)));
% rate_fit=(x(1)*(1/((2*pi)^0.5*x(2)))*exp(-(orien_fit-x(3)).^2./(2*(x(2)).^2))+x(4)*(1/((2*pi)^0.5*x(5)))*exp(-(orien_fit-x(3)-180).^2./(2*(x(5)).^2)));
rate_fit=(x(1)+x(2)*exp(-((orien_fit-x(3)).^2)/(2*x(4)^2))+x(5)*exp(-((orien_fit-x(6)).^2)/(2*x(4)^2)));
% rate_fit=(x(1)+x(2)*exp(-(orien_fit+x(3)).^2./(2*(x(4)).^2)));;

hold on 
plot(orien_fit,rate_fit)
Op=x(3);
