function B0FieldMaps = phase_difference(unwrappedPhase, echoTimes)
%PHASE_DIFFERENCE Computes B0 fieldmaps following the phase difference.
%
% _SYNTAX_
%
%    B0FieldMaps = phase_difference(unwrappedPhase, echoTimes)
%
% _DESCRIPTION_
%
% Returns the B0 field maps that have been calculated from the unwrapped
% phase data corresponding to the different echoes at different acquisition
% times, following the difference phase algorithm.
%
% _INPUT ARGUMENTS_
%
%   unwrappedPhase
%     5D (x,y,z,nEchoes,nAcq) array containing the unwrapped phases.
%
%   echoTime
%     Scalar array (nEchoes, nAcq) containing the different echo times.
%
% _OUTPUTS_
%
%   B0FieldMaps
%     4D (x,y,z,nAcq) array correspondingto the different B0 maps measured
%     at each acquisition time.
%

disp('Computing phase difference B0 maps...')

if size(unwrappedPhase,4) == 1 % If only 1 echo
    echoTimeDiff = echoTimes(1);
    phaseDiff    = unwrappedPhase(:,:,:,:);
else % 2 echoes mapping
    echoTimeDiff = echoTimes(2) - echoTimes(1);
    phaseDiff    = unwrappedPhase(:,:,:,2,:) - unwrappedPhase(:,:,:,1,:);
end

B0FieldMaps = phaseDiff./(2*pi*echoTimeDiff);
