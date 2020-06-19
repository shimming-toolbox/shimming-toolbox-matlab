% b0shim.parts — Shim system components
% 
% Classes in `+b0shim/+parts` define parameter templates for specific shim
% system components. Practically, the classes define properties of the general
% system configurations (namely, config.json files on disk and `b0shim.Config`
% objects in memory). The classes consist only of properties and their
% constructors likely won't need to be called except when creating a new
% configuration.  
% 
% Files
%   Channel - A shim coil element
%   Port    - Serial port protocol parameters
