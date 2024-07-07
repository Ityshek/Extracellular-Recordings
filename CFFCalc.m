function [Data] = CFFCalc(Data,CFFCalcWindow,Raster,SponToResponseFactor,PulseDur,fname,sampling_freq)
%% Initialize Variables
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if answer{1} == '1'
    Data.StimThresh = {}; Data.CFFResponse = []; Data.StimFreq = []; Data.PPC = [];
end
% Find Stim Type
[a,b] = regexp(fname,'_[0123456789]*Hz');
Data.StimFreq = [Data.StimFreq;str2num(fname(a+1:b-2))];
% Calc
CFFResponse = sum(sum(Raster(:,1:CFFCalcWindow)))/size(Raster,1);
CFFSpon = Data.Spon{1}*SponToResponseFactor;
CFFSponStd = Data.SponStd{1}*SponToResponseFactor;
if CFFSpon+2*CFFSponStd < CFFResponse % Check for reponse 2 times greater then std of baseline
    Data.StimThresh = [Data.StimThresh;1];
    Data.CFFResponse = [Data.CFFResponse;CFFResponse];
else
    Data.StimThresh = [Data.StimThresh;0];
    Data.CFFResponse = [Data.CFFResponse;NaN];
end
% PPC Calc
for i=1:64
    PC(i) = PhaseCoherence(i,full(Raster),sampling_freq);
end
Data.PPC = [Data.PPC; PC(str2num(fname(a+1:b-2)))];
figure();plot(linspace(1,64,64),PC);  
xlabel('Frequency [Hz]'); ylabel('PPC Value');
% Plot Final Results
if answer{2} == '1'
    Data.PulseDuration = PulseDur;
    figure();
    plot(Data.StimFreq,Data.PPC)
    
end
ylim([0 1]);
