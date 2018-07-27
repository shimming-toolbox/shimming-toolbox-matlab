This is a WIP for the acdc communication protocol with matlab.

1. Load the arduino code located in matlab_acdc_comm/matlab_acdc_comm.ino. Beaware that it will overwrite your current arduino software.
2. Launch Matlab and call the function matlab_interface_acdc.m
3. Send command on serial port in the form of: a XXX, where XXX is a float. If XXX > 100 the led turns on, else the led turns off
