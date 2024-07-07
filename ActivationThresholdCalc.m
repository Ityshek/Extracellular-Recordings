function [ActiveThresh] = ActivationThresholdCalc(Data,flag)
% Calculate the activation threshold as 2*std of the baseline activity 
if flag == 1 % Check for Alternatig / CFF condition
    for i=1:size(Data.Rast_sort,2)
        for j=1:length(Data.Rast_sort{i})
            Raster = full(Data.Rast_sort{i}{j});
            SponWindow = [35200:43999]; SponWindow2 = [79200:87999]; % Define time windows to search for baseline activity
            for t=1:size(Raster,1)
                CountSpikes1(t) = sum(Raster(t,SponWindow)); CountSpikes2(t) = sum(Raster(t,SponWindow2));
            end
            SponRate(i,j) = mean([CountSpikes2,CountSpikes1])*5; % Calc the mean baseline activiy and normalize to spikes/sec
            SponStd(i,j) = std([CountSpikes2,CountSpikes1])*5; % Calc the std of baseline activiy and normalize to spikes/sec
            ActiveThresh(i,j) = SponRate(i,j)+1*SponStd(i,j); % Define the threshold
        end
    end
else
    for i=1:size(Data.Rast_sort,2)
        Raster = full(Data.Rast_sort{i});
        %SponWindow = [66000:87999];  % Define time windows to search for baseline activity. last 0.5sec in each rep in the CFF case
        SponWindow = [22000:43999];  % Define time windows to search for baseline activity. last 0.5sec in each rep in the FF case
        for t=1:size(Raster,1)
            CountSpikes1(t) = sum(Raster(t,SponWindow)); 
        end
        SponRate(i) = mean(CountSpikes1)*2; % Calc the mean baseline activiy and normalize to spikes/sec
        SponStd(i) = std(CountSpikes1)*2; % Calc the std of baseline activiy and normalize to spikes/sec
        ActiveThresh(i) = SponRate(i)+2*SponStd(i); % Define the threshold
    end
end
end