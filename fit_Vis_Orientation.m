clear all
clc
close all
[filename, pathname]=uigetfile();
openfig([pathname filename ]);
h=get(gca,'children');
orientation=get(h,'xdata');
rate=get(h,'ydata');
spon = rate{1};

 %%
% [Rp,idx] = max(rate{2});
% if idx <= 6
% Rn = rate{2}(idx+6);
% else
% Rn = rate{2}(idx-6);
% end
% Op = orientation{1}(idx);
% B = 0; sigma1 = 10; sigma2 = 10;
% for i=1:length(rate{2})
% R(i) = B+Rp*exp(-((orientation{1}(i)-Op)^2/(2*sigma1^2)))+Rn*exp(-((orientation{1}(i)-Op+180)^2/(2*sigma2^2)));
% end
%% fit the data to a sigmoid  {\displaystyle \Gamma _{\sigma _{1},\sigma _{2}}(x)=I*{\frac {1}{\sigma _{1}{\sqrt {2\pi }}}}\,e^{-{\frac {x^{2}}{2\sigma _{1}^{2}}}}-I*{\frac {1}{\sigma _{2}{\sqrt {2\pi }}}}\,e^{-{\frac {x^{2}}{2\sigma _{2}^{2}}}}.}
figure; 
plot(orientation{1},rate{2},'*k')
hold on

%%
fun=@(x)sum((rate{2}-(x(1)+x(2)*exp(-((orientation{1}-x(3)).^2)/(2*x(4)^2))+x(5)*exp(-((orientation{1}-x(6)).^2)/(2*x(4)^2)))).^2);

 [x eval]=fminsearch(fun, [10 10 60 10 12 60]);

orien_fit=linspace(orientation{1}(1),orientation{1}(end),100);
rate_fit=(x(1)+x(2)*exp(-((orien_fit-x(3)).^2)/(2*x(4)^2))+x(5)*exp(-((orien_fit-x(6)).^2)/(2*x(4)^2)));
% rate_fit=(x(1)+x(2)*exp(-(orien_fit+x(3)).^2./(2*(x(4)).^2)));;

hold on 
plot(orien_fit,rate_fit,'b')
plot(orien_fit,spon(1)*ones(1,length(orien_fit)),'r--')
xlabel('Orientation [Deg]');
xticks([0:60:300]);
xlim([-1 331]);
ylabel('Firing Rate [Spikes/s]');
ylim([0 max(rate_fit)+10])
Op=x(3);
HWHH = round(1.18*x(4),1);