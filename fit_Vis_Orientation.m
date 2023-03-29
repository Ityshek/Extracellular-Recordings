clear all
clc
close all
[filename, pathname]=uigetfile();
openfig([pathname filename ]);
h=get(gca,'children');
orientation=get(h,'xdata');
rate=get(h,'ydata');
spon = rate{1};
count = rate{2} - spon;

for i = 1:length(count)
    if count(i)<0
        count(i) = 0;
    end
end
ScaleCount = count/max(count);
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
close all
figure; 
plot(orientation{1},ScaleCount,'*k')
hold on

%%
fun=@(x)sum((ScaleCount-(x(1)+x(2)*exp(-((orientation{1}-x(3)).^2)/(2*x(4)^2))+x(5)*exp(-((orientation{1}-x(6)).^2)/(2*x(4)^2)))).^2);

 [x eval]=fminsearch(fun, [1 1 10 100 1 0.5]);

orien_fit=linspace(orientation{1}(1),orientation{1}(end),100);
rate_fit=(x(1)+x(2)*exp(-((orien_fit-x(3)).^2)/(2*x(4)^2))+x(5)*exp(-((orien_fit-x(6)).^2)/(2*x(4)^2)));
% rate_fit=(x(1)+x(2)*exp(-(orien_fit+x(3)).^2./(2*(x(4)).^2)));;

hold on 
plot(orien_fit,rate_fit,'b')
xlabel('Orientation [Deg]');
xticks([0:60:300]);
xlim([-1 331]);
ylabel('Scaled Spike Count');
ylim([0 max(rate_fit)])
Op=x(3);
HWHH = round(1.18*x(4),1);