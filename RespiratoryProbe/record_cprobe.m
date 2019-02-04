clear all;
close all;
delete(instrfindall);

%% --
treshold_freq = 100;  % threshold value above/below the mean of previous time points
%% --
duration = 900; % DURATION OF THE MEASUREMENT IN SECONDS

s_cprobe = serial('/dev/cu.usbmodem1538791'); %cprobe 3.5\
s_cprobe.Baudrate = 115200;

%DEFAULT_RESULTS_FOLDER_ROOT = '/Users/alfoi/Desktop/results_cprobe/';
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
figure

FS = stoploop({'Stop cprobe'}) ; 

while(~FS.Stop())
    buffer_data_cprobe = fscanf(s_cprobe,'%s');
    data_cprobe(timepoint) = str2num(buffer_data_cprobe)/100;
    if timepoint ~= 1
        if data_cprobe(timepoint) > data_cprobe(timepoint-1) + treshold_freq || data_cprobe(timepoint) < data_cprobe(timepoint-1) - treshold_freq 
            data_cprobe(timepoint) = data_cprobe(timepoint-1);
        end
    end
    
    time_data(timepoint) = timepoint * 0.1;
    
    plot(time_data,data_cprobe,'r')
    grid on
    title('cprobe')
    xlabel('Time [s]')
    ylabel('Frequency [KHz]')
    drawnow
    
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