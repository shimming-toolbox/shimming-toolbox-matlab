function B0FieldMaps = computeFieldMaps(mappingAlgorithm, unwrappedPhase, phaseJson)
%computeFieldMaps Computes B0 fieldmaps following the specified algorithm
%
%     B0FieldMaps = computeFieldMaps(mappingAlgorithm, unwrappedPhase, phaseJson)
%
% The input `mappingAlgorithm` specifies the way the B0 maps are computed.
% The available algorithms are:    'phaseDifference'
% The input `unwrappedPhase` is the 5D (x,y,z,nEchos,nAcquisitions) output
% of the `unwrap` function and `phaseJson` is the third output of the 
% function `load_niftis` that contains the field `.EchoTime`.
%
% The output B0FielsMaps is a 4D (x,y,z,nAcquisitions) array corresponding
% to the different B0 maps measured at each acquisition time.

switch mappingAlgorithm
    case 'phaseDifference'
        disp('Computing phase difference B0 maps...')
        
        % Different process if only 1 echo
        if size(unwrappedPhase,4) == 1
            echoTimeDiff = phaseJson(1).EchoTime;
            phaseDiff    = unwrappedPhase(:,:,:,:);
        else
            echoTimeDiff = phaseJson(2).EchoTime - phaseJson(1).EchoTime;
            % if using wrapped phase % phaseDiff = angle( wrappedPhase(1) .* conj(wrappedPhase(2) ) ) ; then unwrap
            phaseDiff    = unwrappedPhase(:,:,:,2,:) - unwrappedPhase(:,:,:,1,:);
        end
        
        B0FieldMaps = phaseDiff./(2*pi*echoTimeDiff);
otherwise
    disp(['Unknown mapping algorithm. The available algorithms are:' newline 'phaseDifference'])
end
disp('B0 mapping done')