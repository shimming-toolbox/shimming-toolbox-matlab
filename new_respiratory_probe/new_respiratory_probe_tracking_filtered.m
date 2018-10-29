

clear all;
close all;
delete(instrfindall);

%% --
treshold_freq = 0.3;  % threshold value above/below the mean of previous time points
%% -- 
duration = 200; % DURATION OF THE MEASUREMENT IN SECONDS

s = serial('/dev/cu.usbmodem4471891'); %Teensy 3.5\

exp_descr = inputdlg('Experiment_description: ');
DEFAULT_FREQLOGFILENAME   = ['./results/freqLog-' exp_descr{1} '.bin'] ;
DEFAULT_SAMPLETIMESFILENAME   = ['./results/sampleTimes-' exp_descr{1} '.bin' ] ;
DEFAULT_FREQPLOTFILENAME   = ['./results/freqplot-' exp_descr{1} '.png' ] ;

fopen(s)

timepoint = 1;
limit = duration / 0.1;
figure

while timepoint<limit
    
    data_buffer = fscanf(s,'%s');
    data_freq(timepoint) = str2num(data_buffer)/100;
    % filter the peaks from RF by replacing current value with previous
    % value
    if timepoint ~= 1
        if data_freq(timepoint) > data_freq(timepoint-1) + treshold_freq || data_freq(timepoint) < data_freq(timepoint-1) - treshold_freq 
            data_freq(timepoint) = data_freq(timepoint-1);
        end
    end
    time_data(timepoint) = timepoint * 0.1;
    plot(time_data,data_freq)
    xlabel('Time [s]')
    ylabel('Frequency [kHz]')
    title(strcat('Respiratory_sensor-',exp_descr),'interpreter', 'none');
    grid on
    drawnow
    timepoint = timepoint + 1;
    
end

%% Save freq trace
freqLogFid = fopen( DEFAULT_FREQLOGFILENAME, 'w+' ) ;
fwrite( freqLogFid, data_freq, 'double' ) ;
fclose( freqLogFid );

sampleTimesFid = fopen( DEFAULT_SAMPLETIMESFILENAME, 'w+' ) ;
fwrite( sampleTimesFid, time_data, 'double' ) ;
fclose( sampleTimesFid );

saveas(gcf,DEFAULT_FREQPLOTFILENAME);


fclose(s);