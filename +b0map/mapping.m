function b0FieldMaps = mapping(unwrappedPhase, echoTimes, mappingFunction)
%MAPPING Computes B0 fieldmaps following the specified algorithm.
%
% SYNTAX
%
%    [B0FieldMaps] = mapping(unwrappedPhase, echoTimes)
%    [B0FieldMaps] = mapping(unwrappedPhase, echoTimes, mappingFunction)
%
% DESCRIPTION
%
% Returns the B0 field maps that have been calculated from the unwrapped
% phase data corresponding to the different echoes at different acquisition
% times, following the specified mapping algorithm.
%
% INPUT ARGUMENTS
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
% OUTPUTS
%
%   b0FieldMaps
%     4D (x,y,z,nAcq) array correspondingto the different B0 maps measured
%     at each acquisition time.

narginchk(2,3)

if nargin == 2 % Check the number of arguments
    mappingFunction = 'phase_difference'; % Default mapping function
end

% Check if the specified mapping function exists
if exist(['+b0map/+mappers/' mappingFunction]) ~= 2 
    % If it doesn't, return the list of the available functions
    availableFunctions = string({meta.package.fromName('b0map.mappers').FunctionList.Name});
    error(strjoin(['Mapping function not found. The available functions are:',...
        availableFunctions],'\n'));
end

disp(['Computing the B0 maps using ' mappingFunction '...'])

mapper = str2func(['b0map.mappers.' mappingFunction]);

b0FieldMaps = mapper(unwrappedPhase,  echoTimes);

disp('B0 mapping done')
