function [ PHASE_OUT ] = rescale_phase_image( PHASE_IN )
% rescale siemens phase image (ranging from [0:4095] to the [0:2*PI] range

PHASE_OUT = (PHASE_IN / 4096) * 2*pi;


