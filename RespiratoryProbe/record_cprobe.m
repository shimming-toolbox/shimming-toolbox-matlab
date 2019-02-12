clear all;
close all;
delete(instrfindall);

%% --
treshold_freq = 100;  % threshold value above/below the mean of previous time points
%% --
duration = 900; % DURATION OF THE MEASUREMENT IN SECONDS
port_teensy = ls ('/dev/cu.usbmodem*');

s_cprobe = serial(strcat(port_teensy)) %cprobe 3.5\
%s_cprobe = serial('/dev/tty.usbmodem15387901')
s_cprobe.Baudrate = 115200

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

DEFAULT_CPROBELOGFILENAME   = [DEFAULT_RESULTS_FOLDER '/CProbeLog-' exp_descr{1} '.bin'] ;
DEFAULT_PPROBELOGFILENAME   = [DEFAULT_RESULTS_FOLDER '/PProbeLog-' exp_descr{1} '.bin'] ;
DEFAULT_SAMPLETIMESFILENAME   = [DEFAULT_RESULTS_FOLDER '/sampleTimes-' exp_descr{1} '.bin' ] ;
DEFAULT_PLOTFILENAME   = [DEFAULT_RESULTS_FOLDER '/c_p_probes_plot-' exp_descr{1} '.png' ] ;
DEFAULT_MATFILENAME = [DEFAULT_RESULTS_FOLDER '/variablesDump-' exp_descr{1} '.mat' ] ;

fopen(s_cprobe)

timepoint = 1;

FS = stoploop({'Stop cprobe'}) ;
smooth_window=5;

while(~FS.Stop())
    
    buffer_data_cprobe = fscanf(s_cprobe,'%s');
    time_data(timepoint) = timepoint * 0.1;
    data_cprobe(timepoint) = str2num(buffer_data_cprobe)/100;
    if timepoint ~= 1
        if data_cprobe(timepoint) > data_cprobe(timepoint-1) + treshold_freq || data_cprobe(timepoint) < data_cprobe(timepoint-1) - treshold_freq
            data_cprobe(timepoint) = data_cprobe(timepoint-1);
        end
    end
    
    if timepoint >= 5
        data_cprobe_filtered(timepoint) = mean(data_cprobe((timepoint-smooth_window+1):timepoint))
        figure(1)
        subplot(3,1,1)
        plot(time_data(smooth_window:timepoint),data_cprobe_filtered(smooth_window:timepoint),'r');
        subplot(3,1,2)
        plot(time_data(smooth_window:timepoint),data_cprobe(smooth_window:timepoint),'b');
        grid on
        title('cprobe')
        xlabel('Time [s]')
        ylabel('Frequency [KHz]')
        subplot(3,1,3)
        plot(time_data(smooth_window:timepoint),data_cprobe(smooth_window:timepoint)-data_cprobe_filtered(smooth_window:timepoint),'b');
        drawnow
    end   
    
    timepoint = timepoint + 1;
end

% %% Save freq trace
cprobeLogFid = fopen( DEFAULT_CPROBELOGFILENAME, 'w+' ) ;
fwrite( cprobeLogFid, data_cprobe, 'double' ) ;
fclose( cprobeLogFid );

sampleTimesFid = fopen( DEFAULT_SAMPLETIMESFILENAME, 'w+' ) ;
fwrite( sampleTimesFid, time_data, 'double' ) ;
fclose( sampleTimesFid );

saveas(gcf,DEFAULT_PLOTFILENAME);
save(DEFAULT_MATFILENAME);

fclose(s_cprobe);