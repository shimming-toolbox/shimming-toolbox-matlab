%% IMPORTANT
%Add folders & subfolders to the path the realtime_shimming cloned repo

%%
clc; clear all;close all;

%load('/Users/alfoi/Desktop/work/poly_fit_investigation/results_15-Feb-2019/variablesDump-thomas_base.mat');

%Load frequency measurement file *.mat
uigetfile;

smooth_window=300; % no of samples

figure(1)
timepoint = 1;
FS = stoploop({'Stop cprobe'}) ;
while(~FS.Stop())
    % Wait a number of samples (smooth_window) before fitting a polynom
    if timepoint >smooth_window
        
        %Figure 1  - complete signal filtering
        figure (1)
        
        %Compute polyfit from the beginning of the measurement to timepoint
        subplot(2,1,1)
        p = polyfit(time_data(1:timepoint),data_cprobe(1:timepoint),4);
        y = polyval(p,time_data(1:timepoint));
        plot(time_data(1:timepoint),y,'r');
        
        
        subplot(2,1,2)
        %Substract polyfit from initial signal
        filtered_signal = data_cprobe(1:timepoint) - y;
        plot(time_data(1:timepoint),filtered_signal, 'k');
        title('signal - polyfit 4')
        xlabel('Time [s]')
        ylabel('Frequency [KHz]')
        grid on
        
        %Search for values that pass thru 0 from beginning of the measurement to smooth_window
        if timepoint== smooth_window
            zero_value_locations = find(filtered_signal < 0.01) % this doesn't seem to work
        end
        
    end
    timepoint = timepoint + 1; 
end


