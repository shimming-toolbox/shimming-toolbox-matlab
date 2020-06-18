function B0FieldMaps = field(unwrappedPhase, phaseJson, mappingAlgorithm)
% COMPUTE_FIELD_MAPS Computes B0 fieldmaps following the specified algorithm
%
% _SYNTAX_
%
%    [B0FieldMaps] = compute_Field_Maps(unwrappedPhase, phaseJson)
%    [B0FieldMaps] = compute_Field_Maps(unwrappedPhase, phaseJson, mappingAlgorithm)
%
% _DESCRIPTION_
%
% Returns the B0 field maps that have been calculated from the unwrapped 
% phase data corresponding to the different echoes at different acquisition
% times, following the specified mapping algorithm.
%
% _INPUT ARGUMENTS_
%
%   unwrappedPhase
%     5D (x,y,z,nEchoes,nAcq) array containing the unwrapped phases.
%
%   phaseJson 
%     Structure with the field `.EchoTime` that is an array containing the
%     different echo times corresponding to the measured phases.
%
%   mappingAlgorithm
%     Specifies the algorithm that will be used to compute tht B0 maps.
%     Options are: 'phaseDifference' (default)
%
% _OUTPUTS_
%
%   B0FieldMaps
%     4D (x,y,z,nAcq) array correspondingto the different B0 maps measured
%     at each acquisition time.
%

if nargin == 2 
    mappingAlgorithm = 'phaseDifference';
end

switch mappingAlgorithm
    case 'phaseDifference'
        disp('Computing phase difference B0 maps...')
        
        if size(unwrappedPhase,4) == 1 % Different process if only 1 echo
            echoTimeDiff = phaseJson(1).EchoTime;
            phaseDiff    = unwrappedPhase(:,:,:,:);
        else
            echoTimeDiff = phaseJson(2).EchoTime - phaseJson(1).EchoTime;
            phaseDiff    = unwrappedPhase(:,:,:,2,:) - unwrappedPhase(:,:,:,1,:);
        end
        
        B0FieldMaps = phaseDiff./(2*pi*echoTimeDiff);
        disp('B0 mapping done')
        
otherwise
    disp(['Unknown mapping algorithm. The available algorithms are:' newline '- phaseDifference'])
end
