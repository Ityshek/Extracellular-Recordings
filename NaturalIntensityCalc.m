function [Data] = NaturalIntensityCalc(Data,PSTHbinsize,fname)
%% Natural Intensity Response Curve
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if answer{1} == '1'
    Data.NaturalIntensity = []; Data.NaturalIntensityResponse = cell(1,length(Data.Clusters));   Data.NaturalLatency = cell(1,length(Data.Clusters));
end

% Calculation
% Amplitude to Intensity Conversion
B = [10,20,30,40,50,60,70,80,90,100;0.55,1.2,2.9,5.6,7.9,11.4,18.6,25.6,34.6,42.7]; % 1st row are current in Amp. Second row are intensity in CD/m^2.
a = [strfind(fname,'Hz')+1,strfind(fname,'amp')];
Amp = fname(a(1)+2:a(2)-1);
Data.NaturalIntensity = [Data.NaturalIntensity;B(find(B == str2num(Amp))+1)]; %Convert to Intensity/mm^2 for 40% duty cycle.
% Response Calculation
for i=1:length(Data.Clusters)
    ResponseWindow = [(20/PSTHbinsize)+1:(100/PSTHbinsize)]; % Define time window for Prosthetic response (10-100ms post trigger). add 2 bins for -10 and 0 bins in PSTH.
    Data.NaturalIntensityResponse{i} = [Data.NaturalIntensityResponse{i};round(max(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow)),2)];
end

% Prosthetic Response Latency
for i=1:length(Data.Clusters)
    V = min(find(round(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow),2) == Data.NaturalIntensityResponse{i}(end))+min(ResponseWindow)-1);
    Data.NaturalLatency{i} = [Data.NaturalLatency{i};(V-2)*PSTHbinsize]; % Calculate latency of response in mSec.
end
if answer{2} =='1'
    % Plotting
    %   Spiking Rate
    for i=1:length(Data.Clusters)
        Data.IntensityCurve{i} = figure();
        spon = ones(1,length(Data.NaturalIntensityResponse{i}))*(mean(Data.Spon{i})+(3*std(Data.Spon{i})));
        plot([Data.NaturalIntensity],[Data.NaturalIntensityResponse{i}],'-',Data.NaturalIntensity,spon,'--')
        xticks(linspace(0,round(max(Data.NaturalIntensity),1),5))
        ylabel('Spiking Rate[Hz]','FontSize',20)
        xlabel('Intensity[CD/m^2]','FontSize',20)
        legend('Intensity Response')%,'Spontaneous Activity');
        title(['Unit ',num2str(i)]);
        ylim([max(Data.Spon{i})-5 max(Data.NaturalIntensityResponse{i})+5])
        %   Latency
        Data.LatencyCurve{i} = figure();
        plot(Data.NaturalIntensity,Data.NaturalLatency{i})
        xticks(linspace(0,round(max(Data.NaturalIntensity),1),5))
        ylabel('Latency[ms]','FontSize',20)
        xlabel('Intensity[CD/m^2]','FontSize',20)
        legend('Latency of Response');
        title(['Unit ',num2str(i)]);
        ylim([min(Data.NaturalLatency{i})-10 max(Data.NaturalLatency{i})+10])
    end
end
end
