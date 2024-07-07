% Plot Raster and PSTH of Processed Data
clear all
[file,path]=uigetfile('*.mat','Select a Folder to Load Fig Files From');
load([path,file]);
close all;
figure();
if ~isnan(strfind(file,'Alt'))
    binsize_sec = 0.01; sampling_freq = 44000;
    t_pst=1000*linspace(-10*10^-3,size(Data.Psth_sort{1}{1},2)*binsize_sec-10*10^-3,length(Data.Psth_sort{1}{1}));
    xtickc = linspace(-440,max(t_pst)*44,3);
    names= round(((xtickc/44000)-10^-3),2)*10^3; names(2:end-1) = names(2:end-1) + 10;
    t_Rast = linspace(-440,max(xtickc),size(Data.Rast_sort{1},2));
    NumCPD = length(Data.Rast_sort);
    NumOrs = length(Data.Rast_sort{1});
    t = tiledlayout(NumCPD,NumOrs);
    for i=1:NumCPD
        for j = 1:NumOrs
            flag=1;
            ActiveThresh = ActivationThresholdCalc(Data,flag);
            nexttile
            %subplot(1,NumTrails,i)
            yyaxis left
            spy(Data.Rast_sort{i}{j},6,'k')
            ylabel([num2str(round(Data.CPM(i),0)),'CPM'],'FontSize',10)
            axis square
            %xlim([-440 21560]);
            xtickc=round(linspace(0,round(length(Data.Rast_sort{i}{j})),3),0); names = round(xtickc/44)/1000;
            %    xtickc=round(linspace(-10*10^-3*sampling_freq,round(length(Data.Rast_sort{i})) ...
            %    -10*10^-3*sampling_freq,3),0);
            %   names= round(((xtickc/44)-10^-3),0); names(2:end) = names(2:end)+10;
            set(gca, 'XTick',  xtickc, 'XTickLabel', names)
            set(gca,'YColor','k')
            xlabel([])
            ylim([0 40])
            hold on
            %subplot(1,NumTrails,i)
            yyaxis right
            plot(gca,t_pst*44,Data.Psth_sort{i}{j}*100,'b-','LineWidth',1.5)
            hold on
            plot(gca,t_pst*44,ones(1,length(t_pst))*ActiveThresh(i,j),"Color",'r','LineStyle','-','LineWidth',1)
            % figure settings
            set(gca, 'XTick',  xtickc, 'XTickLabel', names,'FontSize',12,'XTickLabelRotation',0)
            set(gca,'YColor','r')
            xlim([-440 max(xtickc)]);
            %xlim([-440 21560]);
            ylim([0 100]); %ylabel('Spikes/s','FontSize',10);
            %         xlabel('Time [sec]','FontSize',16)
            %         ylabel('Spikes/Sec','FontSize',16)
            title([num2str(Data.Orientations(j)),'^0'],'FontSize',12);
            box off
            set(gca,'color','none')
        end
    end
    %xlabel(t,'Time [sec]','FontSize',16);
    OrsFig = figure(); 
    
    % Build the Combined Response Matrix
    MeanResMat = (Data.AltResponse +Data.Response)/2;
    %[MaxResponse{1},MaxIdx{1}] = max(Data.Response(:)); [MaxResponse{2},MaxIdx{2}] = max(Data.AltResponse(:)); 
    [MaxVal,MaxInd] = max(MeanResMat(:));
    [MaxRow,MaxCol] = ind2sub(size(MeanResMat),MaxInd); 
    
    % Plot all orientations across the CPM
    Ors = rmmissing(split(num2str(Data.Orientations),' '));
    CPMFig = figure();
    for C = 1:NumOrs
    CFig{C} = plot(Data.CPM,MeanResMat(:,C)/MaxVal); 
    hold on
    Ors{C} = [Ors{C},char(176)];
    end
    xlabel('CPM'); ylabel('Normalized Response');
    legend([CFig{1:end}],Ors{:},EdgeColor='none',Color='none');
    axis square; box off;
    set(gca,'color','none','FontSize',15)
    
    % Plot all CPMs across the orientations
    CPMs = rmmissing(split(num2str(round(Data.CPM,0)),' '));
    OrsFig = figure();
    for O = 1:NumCPD
    OFig{O} = plot(Data.Orientations,MeanResMat(O,:)/MaxVal); 
    hold on
    CPMs{O} = [CPMs{O},'CPM'];
    end
    xlabel(['Orientation[',char(176),']']); ylabel('Normalized Response');
    legend([OFig{1:end}],CPMs{:},EdgeColor='none',Color='none');
    axis square; box off;
    set(gca,'color','none','FontSize',15)

else
    binsize_sec = 0.01; sampling_freq = 44000;
    t_pst=1000*linspace(-10*10^-3,size(Data.Psth_sort{1},2)*binsize_sec-10*10^-3,length(Data.Psth_sort{1}));
    xtickc = linspace(-440,max(t_pst)*44,3);
    names= round(((xtickc/44000)-10^-3),2)*10^3; names(2:end-1) = names(2:end-1) + 10;
    t_Rast = linspace(-440,max(xtickc),size(Data.Rast_sort{1},2));
    NumTrails = length(Data.Rast_sort);
    for i=1:NumTrails
        flag=0;
        ActiveThresh = ActivationThresholdCalc(Data,flag);
        nexttile
        %subplot(1,NumTrails,i)
        yyaxis left
        spy(Data.Rast_sort{i},6,'K')
        %ylabel('Repetition','FontSize',12)
        axis square
        xlim([-440 22000]);
        xlabel([])
        xtickc=round(linspace(0,round(length(Data.Rast_sort{i})),5),0); names = round(xtickc/44)/1000;
        %    xtickc=round(linspace(-10*10^-3*sampling_freq,round(length(Data.Rast_sort{i})) ...
        %    -10*10^-3*sampling_freq,3),0);
        %   names= round(((xtickc/44)-10^-3),0); names(2:end) = names(2:end)+10;
        set(gca, 'XTick',  xtickc, 'XTickLabel', names)
        set(gca,'YColor','k')
        ylim([0 40])
        hold on
        %subplot(1,NumTrails,i)
        yyaxis right
        plot(gca,t_pst*44,Data.Psth_sort{i}*100,'b-','LineWidth',1)
        hold on
        plot(gca,t_pst*44,ones(1,length(t_pst))*ActiveThresh(i),"Color",'r','LineStyle','-','LineWidth',1)
        % figure settings
        set(gca, 'XTick',  xtickc, 'XTickLabel', names,'FontSize',12,'XTickLabelRotation',0)
        set(gca,'YColor','r')
        xlim([-440 max(xtickc)]);
        %xlim([-440 44000]);
        ylim([0 100])
        %xlabel('Time [sec]','FontSize',16)
        %ylabel('Spikes/Sec','FontSize',16)
        title(['Pulse Freq: ',num2str(Data.StimFreq(i)),'Hz'],'FontSize',16)
        %title(['Intensity: ',num2str(Data.ProstheticIntensity(i)),'mW/mm^2'],'FontSize',10)
        %title(['Intensity: ',num2str(Data.NaturalIntensity(i)),'nW/mm^2'],'FontSize',10)
        box off
        set(gca,'color','none')
        % figure settings
        %    xlabel('Time [ms]','FontSize',10)
        %     ylabel('Spike Count','FontSize',10)
        %     ylim([0 1]);
        %     xlim([-10 max(t_pst)]);
        %     axis square
        %     hold on
    end
end