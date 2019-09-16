function load_PMU_resp(fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function load_PMU_resp(fname);
%
% Description: A loader for the Respiration data coming off of the Siemens PMU 
%              Based on: https://github.com/timothyv/Physiological-Log-Extraction-for-Modeling--PhLEM--Toolbox
%              
%  
% Input: fname = name of text file with respiration information
%
% Author: Eva Alonso Ortiz, Ecole Polytechnique Montreal
%         eva.alonso.ortiz@gmail.com
%
% Date: Sept. 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[filepath,name,ext] = fileparts(fname);

% Sampling Frequency
Hz = 400;   


fclose('all');fid=fopen(fname);
ignore=textscan(fid,'%s',8); %Ignore first 8 values.

data_blk1 = textscan(fid,'%u16'); %Read data until end of u16 data.

start_data_blk2 = size(data_blk1{1},1);
skip = start_data_blk2+19;

fclose('all');fid=fopen(fname);
ignore = textscan(fid,'%s',skip);

data_blk2 = textscan(fid,'%u16'); %Read data until end of u16 data.

tmp_data_blk1 = data_blk1{1},1;
tmp_data_blk2 = data_blk2{1},1;

data_blk1 = tmp_data_blk1(1:size(tmp_data_blk1)-1);
data_blk2 = tmp_data_blk2(1:size(tmp_data_blk2)-1);

data = vertcat(data_blk1,data_blk2);

footer = textscan(fid,'%s');   %Read in remaining data (time stamps and statistics).


%Get time stamps from footer:
for n=1:size(footer{1},1)
    if strcmp(footer{1}(n),'LogStartMDHTime:')  %log start time
        LogStartTime=str2num(footer{1}{n+1});
    end
    if strcmp(footer{1}(n),'LogStopMDHTime:')   %log stop time
        LogStopTime=str2num(footer{1}{n+1});
    end
    if strcmp(footer{1}(n),'LogStartMPCUTime:') %scan start time
        ScanStartTime=str2num(footer{1}{n+1});
    end
    if strcmp(footer{1}(n),'LogStopMPCUTime:')  %scan stop time
        ScanStopTime=str2num(footer{1}{n+1});
    end
end

% Remove the systems own evaluation of triggers.
t_on  = find(data == 5000);  % System uses identifier 5000 as trigger ON
t_off = find(data == 5003);  % System uses identifier 5003 as trigger OFF


data(t_on) = [];
data(t_off) = [];

signal = data;
time = (1:length(data))./Hz;

LogStartTime_hr = LogStartTime/1000/60/60;
LogStartTime_min = 60*rem(LogStartTime_hr,1);
LogStartTime_sec = 60*rem(LogStartTime_min,1);

figure;
plot(time',signal);
xlabel('Time (s)','FontSize',15);
ylabel('Relative Pressure','FontSize',15);

txt = strcat('LogStartTime is equal to ',num2str(floor(LogStartTime_hr)),' hr ',num2str(floor(LogStartTime_min)),' min ',num2str(floor(LogStartTime_sec)),' sec ');
text(max(time),double(max(signal))+100,txt,'HorizontalAlignment','right','FontSize',15)

fig_name = strcat(name,'_resp_trace');
print('-djpeg',fig_name);

disp(['LogStartTime is equal to ',num2str(floor(LogStartTime_hr)),' hr ',num2str(floor(LogStartTime_min)),' min ',num2str(floor(LogStartTime_sec)),' sec '])

clear;

