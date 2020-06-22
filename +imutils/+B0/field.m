 function B0FieldMaps = field(unwrappedPhase, echoTimes, mappingFunction)
%FIELD Computes B0 fieldmaps following the specified algorithm.
%
% _SYNTAX_
%
%    [B0FieldMaps] = field(unwrappedPhase, echoTimes)
%    [B0FieldMaps] = field(unwrappedPhase, echoTimes, mappingFunction)
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
%   echoTime
%     Scalar array (nEchoes, nAcq) containing the different echo times.
%
%   mappingFunction
%     Specifies the algorithm that will be used to compute tht B0 maps.
%     Options are: 'phase_difference' (default)
%
% _OUTPUTS_
%
%   B0FieldMaps
%     4D (x,y,z,nAcq) array correspondingto the different B0 maps measured
%     at each acquisition time.
%

narginchk(2,3)

if nargin == 2 % Check the number of arguments
    mappingFunction = 'phase_difference'; % Default mapping function
end

% Check if the specified mapping function exists somewhere
if exist(['+imutils/+b0/+mappers/' mappingFunction]) ~= 2 
    % If it doesn't, return the list of the available functions
    error(strjoin(['Mapping function not found. The available functions are:',...
        {meta.package.fromName('imutils.b0.mappers').FunctionList.Name}],'\n'));
end

disp(['Computing B0 maps using' mappingFunction '...'])

mappingFunction = str2func(['imutils.b0.mappers.' mappingFunction]);

B0FieldMaps = mappingFunction(unwrappedPhase,  echoTimes);

disp('B0 mapping done')

