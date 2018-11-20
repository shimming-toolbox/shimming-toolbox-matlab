clear all;
close all;
delete(instrfindall);

%% --
treshold_freq = 1000;  % threshold value above/below the mean of previous time points
%% --
duration = 60; % DURATION OF THE MEASUREMENT IN SECONDS

s_cprobe = serial('/dev/cu.usbmodem4471891'); %cprobe 3.5\
s_cprobe.Baudrate = 115200;

exp_descr = inputdlg('Experiment_description: ');
DEFAULT_CPROBELOGFILENAME   = ['./results/CProbeLog-' exp_descr{1} '.bin'] ;
DEFAULT_PPROBELOGFILENAME   = ['./results/PProbeLog-' exp_descr{1} '.bin'] ;
DEFAULT_SAMPLETIMESFILENAME   = ['./results/sampleTimes-' exp_descr{1} '.bin' ] ;
DEFAULT_PLOTFILENAME   = ['./results/c_p_probes_plot-' exp_descr{1} '.png' ] ;
DEFAULT_MATFILENAME = ['./results/variablesDump-' exp_descr{1} '.mat' ] ;

fopen(s_cprobe)

timepoint = 1;
limit = duration / 0.1;
figure

while timepoint<limit
    fprintf(s_cprobe,'a')
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
    pause(0.01)
    
    
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