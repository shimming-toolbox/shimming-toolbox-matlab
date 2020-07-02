function unwrappedPhase = unwrap_phase(complexArray, unwrapFunction)
%UNWRAP_PHASE Unwraps the phase images following the specified algorithm.
%
% SYNTAX
%
%     [unwrappedPhase] = unwrap_phase(complexArray)
%     [unwrappedPhase] = unwrap_phase(complexArray, unwrapFunction)
%
% DESCRIPTION
%
% Returns the unwrapped phase maps that have been calculated from the
% magnitude and phase data at the different echo times and for each
% acquisition, following the specified unwrapping algorithm.
%
% INPUTs
%
%   complexArray
%     5D (x,y,z,nEchoes,nAcq) complex array containing the B0 magnitude and
%     wrapped phase data.
%
%   unwrapFunction
%     Specifies the algorithm that will be used to compute tht B0 maps.
%     Options are: 'sunwrap' (default)
%
% OUTPUTS
%
%   unwrappedPhase
%     5D (x,y,z,nEchoes,nAcq) array containing the unwrapped phases.

narginchk(1,2)

if nargin == 1
    unwrapFunction = 'sunwrap'; % default unwrapping function
end

% Check if the specified mapping function exists
if exist(['+unwrap/+unwrappers/' unwrapFunction]) ~= 2 
    % If it doesn't, return the list of the available functions
    availableFunctions = string({meta.package.fromName('unwrap.unwrappers').FunctionList.Name});
    error(strjoin(['Unwrapping function not found. The available functions are:',...
        availableFunctions],'\n'));
end

unwrapper= str2func([unwrapFunction]);

disp(['Unwrapping the data using ' unwrapFunction])

for iAcq = 1:size(complexArray,5)
    for iEcho = 1:size(complexArray,4)
        % Unwrap phase using sunwrap
        unwrappedPhase(:,:,:,iEcho,iAcq) = unwrapper(complexArray(:,:,:,iEcho,iAcq), 0.1);
    end
end
        
disp(['Unwrapping done'])
