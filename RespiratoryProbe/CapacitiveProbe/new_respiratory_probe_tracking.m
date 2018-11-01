

clear all;
close all;
delete(instrfindall);

%s = serial('/dev/cu.usbmodem1538791'); %Teensy 3.2
s = serial('/dev/cu.usbmodem4471891'); %Teensy 3.5
duration = 60; % DURATION OF THE MEASUREMENT IN SECONDS

fopen(s)
timepoint = 1;
limit = duration / 0.1;
figure
while timepoint<limit
    data_buffer = fscanf(s,'%s');
    data_freq(timepoint) = str2num(data_buffer)/100;
    time_data(timepoint) = timepoint * 0.1;
    plot(time_data,data_freq)
    xlabel('Time [s]')
    ylabel('Frequency [kHz]')
    title(strcat('Respiratory_sensor_',datestr(now,30)),'interpreter', 'none')
    %axis ([0 duration 50 inf])
    grid on
    drawnow
    %a = a+1;
    timepoint = timepoint + 1;
    
end

DEFAULT_FREQLOGFILENAME   = ['./results/' datestr(now,30) '-freqLog.bin' ] ;
DEFAULT_SAMPLETIMESFILENAME   = ['./results/' datestr(now,30) '-sampleTimes.bin' ] ;
DEFAULT_FREQPLOTFILENAME   = ['./results/' datestr(now,30) '-freqplot.png' ] ;

freqLogFid = fopen( DEFAULT_FREQLOGFILENAME, 'w+' ) ;
fwrite( freqLogFid, data_freq, 'double' ) ;
fclose( freqLogFid );

sampleTimesFid = fopen( DEFAULT_SAMPLETIMESFILENAME, 'w+' ) ;
fwrite( sampleTimesFid, time_data, 'double' ) ;
fclose( sampleTimesFid );

saveas(gcf,DEFAULT_FREQPLOTFILENAME);
fclose(s);