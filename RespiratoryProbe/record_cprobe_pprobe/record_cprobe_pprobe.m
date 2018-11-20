clear all;
close all;
delete(instrfindall);

%% --
treshold_freq = 0.3;  % threshold value above/below the mean of previous time points
%% --
duration = 120; % DURATION OF THE MEASUREMENT IN SECONDS

s_cprobe = serial('/dev/cu.usbmodem4471891'); %cprobe 3.5\
s_pprobe = serial('/dev/tty.usbmodem1421');
s_pprobe.Baudrate = 115200;
s_cprobe.Baudrate = 115200;

exp_descr = inputdlg('Experiment_description: ');
DEFAULT_CPROBELOGFILENAME   = ['./results/CProbeLog-' exp_descr{1} '.bin'] ;
DEFAULT_PPROBELOGFILENAME   = ['./results/PProbeLog-' exp_descr{1} '.bin'] ;
DEFAULT_SAMPLETIMESFILENAME   = ['./results/sampleTimes-' exp_descr{1} '.bin' ] ;
DEFAULT_PLOTFILENAME   = ['./results/c_p_probes_plot-' exp_descr{1} '.png' ] ;
DEFAULT_MATFILENAME = ['./results/variablesDump-' exp_descr{1} '.mat' ] ;

fopen(s_cprobe)
fopen(s_pprobe)


timepoint = 1;
limit = duration / 0.1;
figure

while timepoint<limit
    fprintf(s_cprobe,'a')
    buffer_data_cprobe = fscanf(s_cprobe,'%s');
    data_cprobe(timepoint) = str2num(buffer_data_cprobe);
    
    fprintf(s_pprobe,'a')
    buffer_data_pprobe = fscanf(s_pprobe,'%f');
    data_pprobe(timepoint) = buffer_data_pprobe;
    
    time_data(timepoint) = timepoint * 0.1;
    
    subplot(2,1,1)
    plot(time_data,data_cprobe,'r')
    grid on
    title('cprobe')
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    drawnow
    
    subplot (2,1,2)
    plot(time_data,data_pprobe,'b')
    grid on
    title('pprobe')
    xlabel('Time [s]')
    ylabel('Pressure [a.u.]')
    drawnow
    
    timepoint = timepoint + 1;
    pause(0.01)
    
    
end

% %% Save freq trace
cprobeLogFid = fopen( DEFAULT_CPROBELOGFILENAME, 'w+' ) ;
fwrite( cprobeLogFid, data_cprobe, 'double' ) ;
fclose( cprobeLogFid );

pprobeLogFid = fopen( DEFAULT_PPROBELOGFILENAME, 'w+' ) ;
fwrite( cprobeLogFid, data_pprobe, 'double' ) ;
fclose( cprobeLogFid );

sampleTimesFid = fopen( DEFAULT_SAMPLETIMESFILENAME, 'w+' ) ;
fwrite( sampleTimesFid, time_data, 'double' ) ;
fclose( sampleTimesFid );

saveas(gcf,DEFAULT_PLOTFILENAME);
save(DEFAULT_MATFILENAME);

fclose(s_pprobe);
fclose(s_cprobe);