function [stimulus_times, stimulus_indexes]=find_stim(stim_Data,stim_sampling_rate,sampling_freq,Begin_record)

            %converting trigger indecies to times in secs
%             stimulus_times=(stim_Data(1,temp))/(stim_sampling_rate);
            stimulus_times=((stim_Data(1,1:2:end))/stim_sampling_rate)-Begin_record;
            stimulus_indexes=round((stimulus_times)*sampling_freq);
         

           