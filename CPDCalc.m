function [Data] = CPDCalc(Data,PSTHbinsize,fname,sf)
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if answer{1} == '1'
    Data.CPDs = []; Data.CPDResponse = cell(1,length(Data.Clusters)); Data.CPDCount = cell(1,length(Data.Clusters));
end

%% Calculation
binsize_sec = PSTHbinsize*10^-3;
%ResponseWindow = [1/binsize_sec:1.1/binsize_sec+20]; % Select for max bin calc
ResponseWindow= [1.05/binsize_sec+1:1.2/binsize_sec+1]; % Define 100ms response window for spike counts and rate.
% a = [max(strfind(fname,'0_')),strfind(fname,'CPD')];
% CPD = fname(a(1):a(2)-1);
% CPD(strfind(CPD,'_')) = '.';
[Int1,Int2] = regexp(fname,'_0.\d*CPD');
 Data.CPDs = [Data.CPDs; str2num(fname(Int1+1:Int2-3))];
 for i=1:length(Data.Clusters)
    Data.CPDResponse{i} = [Data.CPDResponse{i}; max(max(Data.Psth_sort{Data.Clusters(i)}(ResponseWindow)))]; % Calculate spike count response window post pattern onset
    CountCalc = sum(sum(Data.Rast_sort{i}(:,min(ResponseWindow)*10^-2*sf:max(ResponseWindow)*10^-2*sf)))...
        /size(Data.Rast_sort{i},1);
    Data.CPDCount{i} = [Data.CPDCount{i};CountCalc];
end

%% Plotting
if answer{2} =='1'
    for i=1:length(Data.Clusters)
        figure();
        spon = ones(1,length(Data.CPDResponse{i}))*mean(Data.Spon{i},'omitnan');     
        plot(Data.CPDs,Data.CPDCount{i},'-',Data.CPDs,Data.Spon{i},'r-');
        hold on
        plot(Data.CPDs,Data.CPDCount{i}-Data.Spon{i})
        xticks(round(linspace(min(Data.CPDs),max(Data.CPDs),10),2));
        ylabel('Spike Counts','FontSize',20);
        xlabel('CPD','FontSize',20);
        legend('Response','Spontaneous Activity','')
        %ylim([0 100]);
    end
end
end