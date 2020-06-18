function [ PHASE_OUT ] = rescalePhaseImage( PHASE_IN, invert )
% rescale phase image outputed by scanner (ranging from [-4096:4095]
% to the [-PI:PI] range
%
% PHASE_IN: phase image with value ranging from [-4096:4095]
% invert: Invert sign (i.e. to convert from left-handed system to
%         right-handed system (default = 1)

if (~exist('invert', 'var'))
    invert = 1;
end

PHASE_OUT = (PHASE_IN / 4096) * pi;

if (invert ~= 0)
    PHASE_OUT = -PHASE_OUT;
end

end

