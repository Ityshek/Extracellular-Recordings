function [Data] = ProstheticIntensityCalc(Data,PSTHbinsize,fname)
%% Prosthetic Intensity Response Curve
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if answer{1} == '1'
    Data.ProstheticIntensity = []; Data.ProstheticIntensityResponse = cell(1,length(Data.Clusters));   Data.ProstheticLatency = cell(1,length(Data.Clusters));
end

% Calculation
% Amplitude to Intensity Conversion
B = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7;0.42*10^(-3),0.3,0.9,1.45,1.93,2.38,2.8,3.15,3.6,4.05,4.47,4.88,5.3,5.83]; % 1st row are current in Amp. Second row are intensity in mW.
D = [0.8,0.95,1.25,1.6,2.1,2.7,3.9,4.9;0.097,0.27,0.53,1.01,1.47,2.09,3.03,4.05];
a = [strfind(fname,'Hz')+1,strfind(fname,'amp')];
Amp = fname(a(1)+2:a(2)-1);
%Data.ProstheticIntensity = [Data.ProstheticIntensity;D(find(D == str2num(Amp))+1)]; %Convert to Intensity/mm^2 for 40% duty cycle.
Data.ProstheticIntensity = [Data.ProstheticIntensity;B(find(B == str2num(Amp))+1)]; %Convert to Intensity/mm^2 for 40% duty cycle.
% Response Calculation
for i=1:length(Data.Clusters)
    ResponseWindow = [(20/PSTHbinsize)+1:(100/PSTHbinsize)]; % Define time window for Prosthetic response (10-100ms post trigger). add 2 bins for -10 and 0 bins in PSTH.
    Data.ProstheticIntensityResponse{i} = [Data.ProstheticIntensityResponse{i};round(max(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow)),2)];
end

% Prosthetic Response Latency
for i=1:length(Data.Clusters)
    V = min(find(round(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow),2) == Data.ProstheticIntensityResponse{i}(end))+min(ResponseWindow)-1);
    Data.ProstheticLatency{i} = [Data.ProstheticLatency{i};(V-2)*PSTHbinsize]; % Calculate latency of response in mSec.
end
if answer{2} =='1'
    % Plotting
    %   Spiking Rate
    for i=1:length(Data.Clusters)
        Data.IntensityCurve{i} = figure();
        spon = ones(1,length(Data.ProstheticIntensityResponse{i}))*(mean(Data.Spon{i})+(3*std(Data.Spon{i})));
        plot([Data.ProstheticIntensity],[Data.ProstheticIntensityResponse{i}],'-',Data.ProstheticIntensity,spon,'--')
        xticks(linspace(0,round(max(Data.ProstheticIntensity),1),5))
        ylabel('Spiking Rate[Hz]','FontSize',20)
        xlabel('Intensity[mW/mm^2]','FontSize',20)
        legend('Intensity Response')%,'Spontaneous Activity');
        title(['Unit ',num2str(i)]);
        ylim([max(Data.Spon{i})-5 max(Data.ProstheticIntensityResponse{i})+5])
        %   Latency
        Data.LatencyCurve{i} = figure();
        plot(Data.ProstheticIntensity,Data.ProstheticLatency{i})
        xticks(linspace(0,round(max(Data.ProstheticIntensity),1),5))
        ylabel('Latency[ms]','FontSize',20)
        xlabel('Intensity[mW/mm^2]','FontSize',20)
        legend('Latency of Response');
        title(['Unit ',num2str(i)]);
        ylim([min(Data.ProstheticLatency{i})-10 max(Data.ProstheticLatency{i})+10])
    end
end
end