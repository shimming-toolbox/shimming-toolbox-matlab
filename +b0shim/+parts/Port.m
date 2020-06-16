classdef Port 
%b0shim.parts.Port Serial port protocol parameters
%    
%    self = b0shim.parts.Port()
% 
% `Port` properties are simply 6 of the 7 the parameters passed to MATLAB's
% `serialport` function when opening the communication line between a PC and
% the microcontroller directing the shim amplifiers. Hence, the properties are
% only relevant for hardware systems permitting direct control via MATLAB—USB
% (i.e. systems that subclass `b0shim.Com` to define a communication interface).
% For a host MRI, or a virtual shim array, this classdef is likely irrelevant.
%
% __NOTES__  
%
% 1. While `Port` properties all correspond to `serialport()` arguments,
% their default values differ (`Port` defaults are informed by the
% microcontrollers in use at our lab).
%
% 2. As a safety feature, the mandatory first argument to `serialport()`—the
% device ID—is not included as a `Port` property: when connecting to a device,
% a user still needs to specify the port name. (The parameter, in any case,
% tends to vary as device names are determined by the OS.)
% 
% See also  
% [ serialport ](https://www.mathworks.com/help/matlab/ref/serialport.html)  
% [ b0shim.parts.Contents ](./Contents.md)  
% [ b0shim.Specs ](../Specs.md)  

properties

    % 2nd argument to serialport
    BaudRate(1,1) uint32 {mustBeInteger,mustBePositive} = 9600;

    % Name-value argument to serialport
    DataBits(1,1) uint8 {mustBeMember(DataBits, [8 7 6 5] )} = 8;

    % Name-value argument to serialport
    StopBits(1,1) single {mustBeMember(StopBits, [1 1.5 2] )} = 1;

    % Name-value argument to serialport
    FlowControl(1,1) string {mustBeMember(FlowControl, ["none" "hardware" "software"] )} = "none";

    % Name-value argument to serialport
    ByteOrder(1,1) string {mustBeMember(ByteOrder, ["little-endian" "big-endian"] )} = "big-endian";

    % Name-value argument to serialport
    Timeout(1,1) single {mustBeNumeric, mustBePositive} = 10;
    
end

% =========================================================================
% =========================================================================
methods 
% =========================================================================
function self = Port( )
    return
end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
