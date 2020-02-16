function [DEFAULTS] = makedefaults( varargin )
%MAKEDEFAULTS Return default parameters for MrdiIo.make()
%
% DEFAULTS     = MrdiIo.makedefaults( )
% defaultValue = MrdiIo.makedefaults( fieldName )
% 
% Called without arguments, the complete default parameters struct is returned.
% If a valid fieldName [char or string] is provided as a single input argument
% then the corresponding default value is returned.
% Valid field names correspond to MrdiIo properties (e.g. 'isSearchRecursive', 'nBytesMax') 
% 
% For explanations of their respective significance, see doc MrdiIo
    
narginchk(0,1) ;

Mc = metaclass( MrdiIo ) ;

for iProp = 1 : length( Mc.PropertyList ) 
    if ~Mc.PropertyList(iProp).Constant && Mc.PropertyList(iProp).HasDefault 
        DEFAULTS.( Mc.PropertyList(iProp).Name ) = Mc.PropertyList(iProp).DefaultValue ;
    end
end

if nargin == 0
    return
elseif nargin == 1 
    validateattributes( varargin{1}, {'string', 'char'}, {'scalartext'}, mfilename, 'searchDir', 1 ) ;
    if isfield( DEFAULTS, varargin{1} )
        DEFAULTS = DEFAULTS.( varargin{1} ) ;
        return ;
    else
        warning('Invalid field name, returning complete default parameters struct.') ; 
        return ;
    end
else
    error('Invalid input') ; 
end

end
