function [MmInfo] = metamethodinfo( Mm )
%METAMETHODINFO Return struct of info derived from a meta.method object array
narginchk(1,1) ;
assert( strcmp( class( Mm ), 'meta.method' ), 'Invalid inputs' ) ;

PROPS = [ "Access" ;
          "Static" ;
          "Abstract" ;
          "Sealed" ;
          "ExplicitConversion" ;
          "Hidden" ;
          "InputNames" ;
          "OutputNames" ;
          "DefiningClass" ;
        ] ;

nObj = numel( Mm ) ;

for iObj = 1 : nObj 

    PROPERTIES_TO_COPY( iObj )
    MmInfo( iObj ).( =

           
           ] ;

    MmInfo = struct ( M

    Name
    % Description
    % DetailedDescription
    Access
    Static
    Abstract
    Sealed
    ExplicitConversion
    Hidden
    InputNames
    OutputNames
    
    
    Mthds = [] ;




% Copies certain elements of a meta.method object's properties to the fields of
% a returned struct array. 
%
% ### Syntax ###
%
%  [MmInfo] = METAMETHODINFO( Mm )

% ### Inputs ###
% 
% - Mc: The meta.class intance. 
%
% - methodNames: string or cellstr of the specific method names of interest.
%   [default: all methods listed in Mc.MethodList)
%
% ### Outputs ###
%
% - Mthds: Struct array with field values copied from Mc.MethodList: 
%
% Name
% Description
% DetailedDescription
% Access
% Static
% Abstract
% Sealed
% ExplicitConversion
% Hidden
% InputNames
% OutputNames
% DefiningClass (class names as strings) 
% 
%  [Mthds] = GETMETHODATTRIBUTES( Mc, methodNames )
%
% ### Notes ### 
%
% 1. Currently, all the hidden properties (and some of the non-hidden
% ones) of meta.method objects appear to be quite useless--despite their
% promising sounding names! Therefore, these are not copied to Att. 
%
% 2. For now, to (hopefully) avoid any strange bugs/behaviour from the
% meta.class package, the list of properties to copy from the meta.method
% instance is explictly/hard-coded rather than copying all the properties
% wholesale.  
%


    if nargin == 2

        mustBeStringOrCharOrCellstr( methodNames ) ;
        methodNames = strip( string( methodNames ) ) ;

        McM = McM( string({McM.Name}) == methodNames ) ;
        assert( ~isempty( McM ), 'Input meta.class instance does not contain any of the specified method names' ) ;
    
    end

    % for McM.DefiningClass, just keep the name:
    Mthds.DefiningClass = string( McM.DefiningClass.Name ) ;
    
    for iField = 1 : numel( FIELDS )
        Mthds.( FIELDS(iField) ) = McM( iField ) ;
    end
    
    Mthds.nInputs  = numel( Mthds.InputNames ) ;    
    Mthds.nOutputs = numel( Mthds.OutputNames ) ;    

end

%GETMETHODATTRIBUTES Return struct array containing methods attributes copied from a meta.method object

% METHODSFROMMETACLASS 
    

