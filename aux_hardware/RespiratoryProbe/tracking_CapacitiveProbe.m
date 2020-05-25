clear all;
close all;
delete(instrfindall);

%% Initalize communication and create folder structure
port_teensy = ls ('/dev/cu.usbmodem*');

s_cprobe = serial(strcat(port_teensy));
s_cprobe.Baudrate = 115200;

DEFAULT_RESULTS_FOLDER_ROOT = '~/Desktop/results_cprobe/';

str_date = date;
DEFAULT_RESULTS_FOLDER = strcat('~/Desktop/results_cprobe/results_',str_date);

if ~exist(DEFAULT_RESULTS_FOLDER_ROOT, 'dir')
    mkdir(DEFAULT_RESULTS_FOLDER_ROOT)
end

if ~exist(DEFAULT_RESULTS_FOLDER, 'dir')
    mkdir(DEFAULT_RESULTS_FOLDER)
end
exp_descr = inputdlg('Experiment_description: ');

dt=datestr(now,'yyyymmdd_HH-MM-SSAM');

DEFAULT_CPROBELOGFILENAME   = [DEFAULT_RESULTS_FOLDER '/' dt '_CProbeLog_' exp_descr{1} '.bin'] ;
DEFAULT_PPROBELOGFILENAME   = [DEFAULT_RESULTS_FOLDER '/' dt '_PProbeLog_' exp_descr{1} '.bin'] ;
DEFAULT_SAMPLETIMESFILENAME   = [DEFAULT_RESULTS_FOLDER '/' dt '_SampleTimes_' exp_descr{1} '.bin' ] ;
DEFAULT_PLOTFILENAME   = [DEFAULT_RESULTS_FOLDER '/' dt '_CPprobesPlot_' exp_descr{1} '.png' ] ;
DEFAULT_MATFILENAME = [DEFAULT_RESULTS_FOLDER '/' dt '_VariablesDump_' exp_descr{1} '.mat' ] ;

fopen(s_cprobe)

timepoint = 1;

FS = stoploop({'Stop cprobe'}) ;
%% No of measurements before the actual polyfit will begin
smooth_window=5;

%% Data acquisition until stop button is pushed
figure(1)
while(~FS.Stop())
    
    buffer_data_cprobe = fscanf(s_cprobe,'%s');
    time_data(timepoint) = timepoint * 0.1;
    data_cprobe(timepoint) = str2num(buffer_data_cprobe)/100;
   
    if timepoint >= smooth_window
        
        subplot(2,1,1)
        % Polyfit 4th degree to get rid of DC
        p = polyfit(time_data(smooth_window:timepoint),data_cprobe(smooth_window:timepoint),4);
        y = polyval(p,time_data(smooth_window:timepoint));
        plot(time_data(smooth_window:timepoint),y,'r');
        hold on
        title('Legend: normal signal(blue), polyfit 4 (red)');
        plot(time_data(smooth_window:timepoint),data_cprobe(smooth_window:timepoint),'b');
        xlabel('Time [s]')
        ylabel('Frequency [KHz]')
        grid on
        hold off
        drawnow
        
        subplot(2,1,2)
        % Filtered signal (signal - polyfit 4)
        data_cprobe_filtered = data_cprobe(smooth_window:timepoint) - y;
        plot(time_data(smooth_window:timepoint), data_cprobe_filtered, 'k');
        title('signal - polyfit 4')
        xlabel('Time [s]')
        ylabel('Frequency [KHz]')
        grid on 
        
    end
    
    timepoint = timepoint + 1;
end

%% Save freq trace
cprobeLogFid = fopen( DEFAULT_CPROBELOGFILENAME, 'w+' ) ;
fwrite( cprobeLogFid, data_cprobe, 'double' ) ;
fclose( cprobeLogFid );

sampleTimesFid = fopen( DEFAULT_SAMPLETIMESFILENAME, 'w+' ) ;
fwrite( sampleTimesFid, time_data, 'double' ) ;
fclose( sampleTimesFid );

saveas(gcf,DEFAULT_PLOTFILENAME);
save(DEFAULT_MATFILENAME);

fclose(s_cprobe);