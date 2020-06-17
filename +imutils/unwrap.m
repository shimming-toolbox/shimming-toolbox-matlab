function unwrappedPhase = unwrap(unwrapAlgorithm, mag, phase)
%unwrap Computes the unwrapped phase following the specified algorithm
%
%     unwrappedPhase = unwrap(unwrapAlgorithm, mag, phase)
%
% The input `unwrapAlgorithm` specifies the way the B0 maps are unwrapped.
% The available algorithms are:    'sunwrap'
% The inputs `mag` ans `phase` are 5D (x,y,z,nEchos,nAcquisitions) arrays.
%
% The output `unwrappedPhase` is a 5D (x,y,z,nEchos,nAcquisitions) array
% containing the results of the unwrapping process.

switch unwrapAlgorithm
    case 'sunwrap'  % If the chosen unwrapper is sunwrap
        disp('Unwrapping with sunwrap...')
        
        % Init Unwrapped Phase
        for iAcq = 1:size(phase,5)
            for iEcho = 1:size(phase,4)
                % Get the magnitude for a specific echo
                magNorm = mat2gray(mag(:,:,:,iEcho,iAcq));
                
                % Calculate the phase in radians, assumes there are wraps
                phasePi = mat2gray(phase(:,:,:,iEcho,iAcq))*2*pi - pi;
                
                % Unwrap phase using sunwrap
                unwrappedPhase(:,:,:,iEcho,iAcq) = sunwrap(magNorm .* exp( 1i* phasePi ), 0.1);
                
            end
        end
    otherwise
        disp(['Unknown algorithm. The available algorithms are:' newline '- sunwrap'])
end
disp('Unwrapping done')