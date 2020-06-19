function unwrappedPhase = unwrap_phase(complexArray, unwrapFunction)
%UNWRAP_PHASE Unwraps the phase images following the specified algorithm.
%
% _SYNTAX_
%
%     [unwrappedPhase] = unwrap_phase(complexArray)
%     [unwrappedPhase] = unwrap_phase(complexArray, unwrapFunction)
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
%   unwrapFunction
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
    unwrapFunction = 'sunwrap'; % default unwrapping function
end
 
% Check if the specified mapping function exists somewhere
if exist(['+imutils/+unwrappers/' unwrapFunction]) ~= 2 
    % If it doesn't, return the list of the available functions
    error(strjoin(['Mapping function not found. The available functions are:',...
        {meta.package.fromName('imutils.unwrappers').FunctionList.Name}],'\n'));
end

mappingFunction = str2func(['imutils.unwrappers.' unwrapFunction]);

unwrappedPhase = mappingFunction(complexArray);
        
disp(['Unwrapping done'])
