function [Data] = CPDCalc(Data,PSTHbinsize,fname,sf,Rast_sort,Psth_sort)
prompt = {'Is First Trail? (1=YES 0=NO)','Is Final Trail? (1=YES 0=NO)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if answer{1} == '1'
    Data.CPDs = []; Data.CPDResponse = []; Data.CPDCount = [];
end

%% Calculation
binsize_sec = PSTHbinsize*10^-3;
%ResponseWindow = [1/binsize_sec:1.1/binsize_sec+20]; % Select for max bin calc
ResponseWindow= [1.05/binsize_sec+1:1.25/binsize_sec+1]; % Define 200ms response window for spike counts and rate.
% a = [max(strfind(fname,'0_')),strfind(fname,'CPD')];
% CPD = fname(a(1):a(2)-1);
% CPD(strfind(CPD,'_')) = '.';
[Int1,Int2] = regexp(fname,'_[.0123456789]*CPD');
 Data.CPDs = [Data.CPDs; str2num(fname(Int1+1:Int2-3))];
Data.CPDResponse = [Data.CPDResponse; max(max(Psth_sort(ResponseWindow)))]; % Calculate spike count response window post pattern onset
CountCalc = sum(sum(Rast_sort(:,min(ResponseWindow)*10^-2*sf:max(ResponseWindow)*10^-2*sf)))...
    /size(Rast_sort,1);
Data.CPDCount = [Data.CPDCount;CountCalc];



%% Plotting
if answer{2} =='1'
        for i=1:length(Data.Rast_sort)
            SponCalc(i) = sum(sum(Data.Rast_sort{i}(:,end-size(ResponseWindow,2)*10^-2*sf:end)))...
    /size(Data.Rast_sort{i},1);
        end
        figure();
        %spon = ones(1,length(Data.CPDResponse{i}))*mean(Data.Spon{i},'omitnan');     
        %plot(Data.CPDs,Data.CPDCount{i},'-',Data.CPDs,Data.Spon{i},'r-');
        %hold on
        CorrectedCalc = Data.CPDCount-SponCalc(i);
        CorrectedCalc(find(CorrectedCalc<0)) = 0;
        plot(Data.CPDs,CorrectedCalc)
        xticks(round(linspace(min(Data.CPDs),max(Data.CPDs),5),1));
        ylabel('Spike Counts','FontSize',20);
        xlabel('CPD','FontSize',20);
        %legend('Response')
        %ylim([0 100]);
ax.PlotBoxAspectRatio = [1,1,1]; ax.FontSize = 20;
ax.Box = 'off'; ax.Color = "none";
end
end