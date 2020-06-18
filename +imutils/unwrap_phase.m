function unwrappedPhase = unwrap_phase(complexArray, unwrapAlgorithm)
%UNWRAP_PHASE Unwraps the phase images following the specified algorithm.
%
% _SYNTAX_
%
%     [unwrappedPhase] = unwrap_phase(complexArray)
%     [unwrappedPhase] = unwrap_phase(complexArray, unwrapAlgorithm)
%
% _DESCRIPTION_
%
% Returns the unwrapped phase maps that have been calculated from the
% magnitude and phase data at the different echo times and for each
% acquisition, following the specified unwrapping algorithm.
%
% _INPUT ARGUMENTS_
%
%   complexArray
%     5D (x,y,z,nEchoes,nAcq) complex array containing the B0 magnitude and
%     wrapped phase data.
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

narginchk(1,2)

if nargin == 1
    unwrapAlgorithm = 'sunwrap';
end

switch unwrapAlgorithm
    case 'sunwrap' % If the chosen unwrapper is sunwrap
        disp('Unwrapping with sunwrap...')
        
        % Init Unwrapped Phase
        for iAcq = 1:size(complexArray,5)
            for iEcho = 1:size(complexArray,4)
                
                % Unwrap phase using sunwrap
                unwrappedPhase(:,:,:,iEcho,iAcq) = sunwrap(complexArray(:,:,:,iEcho,iAcq), 0.1);
                
            end
        end
        disp('Unwrapping done')
        
    otherwise % If unwrapAlgorithm doesn't match any available algorithm
        disp(['Unknown algorithm. The available algorithms are:' ...
            newline '- sunwrap'])
end
