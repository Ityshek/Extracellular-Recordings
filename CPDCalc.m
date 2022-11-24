function [Data] = CPDCalc(Data,PSTHbinsize,fname)
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if answer{1} == '1'
    Data.CPDs = cell(1,length(Data.Clusters)); Data.CPDResponse = cell(1,length(Data.Clusters));
end

%% Calculation

binsize_sec = PSTHbinsize*10^-3;
ResponseWindow = [1/binsize_sec:1.1/binsize_sec+20]; % Select for max bin calc
ResponseWindow= [1/binsize_sec:1.1/binsize_sec+20];
a = [max(strfind(fname,'0_')),strfind(fname,'CPD')];
CPD = fname(a(1):a(2)-1);
CPD(strfind(CPD,'_')) = '.';
for i=1:length(Data.Clusters)
    Data.CPDs{i} = [Data.CPDs{i}; str2num(CPD)];
    Data.CPDResponse{i} = [Data.CPDResponse{i}; max(max(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow)))];
end

%% Plotting
if answer{2} =='1'
    for i=1:length(Data.Clusters)
        figure();
        spon = ones(1,length(Data.CPDResponse{i}))*mean(Data.Spon{i},'omitnan');
        plot(Data.CPDs{i},Data.CPDResponse{i},'-',Data.CPDs{i},spon,'r-');
        xticks(round(linspace(min(Data.CPDs{i}),max(Data.CPDs{i}),10),2));
        ylabel('Firing Rate [Spikes/s]','FontSize',20);
        xlabel('CPD','FontSize',20);
        legend('CPD Response','Spontaneous Activity')
        ylim([0 100]);
    end
end
end