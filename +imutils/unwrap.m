function unwrappedPhase = unwrap(mag, phase, unwrapAlgorithm)
%UNWRAP Unwraps the phase images following the specified algorithm
%
% _SYNTAX_
%
%     [unwrappedPhase] = unwrap(mag, phase)
%     [unwrappedPhase] = unwrap(mag, phase, unwrapAlgorithm)
%
% _DESCRIPTION_
%
% Returns the unwrapped phase maps that have been calculated from the
% magnitude and phase data at the different echo times and for each
% acquisition, following the specified unwrapping algorithm.
%
% _INPUT ARGUMENTS_
%
%   mag
%     5D (x,y,z,nEchoes,nAcq) array containing the magnitude data.
%
%   phase 
%     5D (x,y,z,nEchoes,nAcq) array containing the phase data.
%
%   unwrapAlgorithm
%     Specifies the algorithm that will be used to compute tht B0 maps.
%     Options are: 'sunwrap' (default)
%
% _OUTPUTS_
%
%   B0FieldMaps
%     5D (x,y,z,nEchoes,nAcq) array containing the unwrapped phases.
%

if nargin == 2 
    unwrapAlgorithm = 'sunwrap';
end

switch unwrapAlgorithm
    case 'sunwrap' % If the chosen unwrapper is sunwrap
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
        disp('Unwrapping done')
        
    otherwise
        disp(['Unknown algorithm. The available algorithms are:' newline '- sunwrap'])
end
