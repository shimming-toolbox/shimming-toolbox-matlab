function unwrappedPhase = sunwrap(complexArray)
%SUNWRAP Calls the sunwrap algorithm to unwrap the phase images.
%
% _SYNTAX_
%
%    unwrappedPhase = sunwrap(complexArray)
%
% _DESCRIPTION_
%
% Returns the unwrapped phase maps that have been calculated from the
% magnitude and phase data at the different echo times and for each
% acquisition, using the sunwrap Matlab function.
%
% _INPUT ARGUMENTS_
%
%   complexArray
%     5D (x,y,z,nEchoes,nAcq) complex array containing the B0 magnitude and
%     wrapped phase data.
%
% _OUTPUTS_
%
%   B0FieldMaps
%     5D (x,y,z,nEchoes,nAcq) array containing the unwrapped phases.
%

disp('Unwrapping with sunwrap...')

% Init Unwrapped Phase
for iAcq = 1:size(complexArray,5)
    for iEcho = 1:size(complexArray,4)
        % Unwrap phase using sunwrap
        unwrappedPhase(:,:,:,iEcho,iAcq) = sunwrap(complexArray(:,:,:,iEcho,iAcq), 0.1);
    end
end
