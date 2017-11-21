
% Initialization
clear all;
delete(instrfindall);
% Define commands
acdc_cmds= {'a','b','c','d','e', 'f','g','h','z'};
current_val =[120,100,100,100,100,100,100,100];



%% Initialize Communication port
%==========================================================================

%function [ComPort] = initializecomport()
if ismac
    portName = '/dev/tty.usbmodem1411'; % Alex's 2015 MacBook Pro -- left USB port
    ComPort = serial( portName,'BaudRate', 115200) ;
end
%end

%% Set current
%==========================================================================
%function [] = set_current(current_val)
fopen (ComPort);
pause (2); % pause required to establish communication
%for ii = 1:1:8
ii = 1;
cmd =strcat(acdc_cmds{ii},{' '},num2str(current_val(ii)));
pause(0.1);
disp(cmd);
fprintf(ComPort,'%s',cmd{1},'sync');
a = fscanf(ComPort,'%s');
disp(strcat('Feeedback from serial: ',a));
%end
%end