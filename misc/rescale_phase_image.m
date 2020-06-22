function [ phaseOut ] = rescale_phase_image( phaseIn )
% rescale siemens phase image (ranging from [0:4095] to the [0:2*PI] range

phaseOut = (phaseIn / 4096) * 2*pi;


