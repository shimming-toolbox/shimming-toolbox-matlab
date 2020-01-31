function isSame = comparegrids( varargin )
%COMPAREGRIDS  Return TRUE if voxel positions coincide 
%   
% Wrapper function to isequal().
%
% isSame = COMPAREGRIDS( Grid1, Grid2 )
% isSame = COMPAREGRIDS( xyz1, xyz2 )
% isSame = COMPAREGRIDS( X1, Y1, Z1, X2, Y2, Z2 )
%
% See also isequal

    if nargin == 2
        if isa( varargin{1}, 'MrdiGrid' ) && isa( varargin{2}, 'MrdiGrid' )
            isSame = ( varargin{1} == varargin{2} ) ;   
        elseif isnumeric( varargin{1} ) && isnumeric( varargin{2} )
            isSame = isequal( varargin{1}, varargin{2} ) ; 
        end
    elseif ( nargin ==6 ) && all( cellfun( @isnumeric, varargin ) )
        isSame  = isequal( varargin{1}, varargin{4} ) ... % compare X-coordinates
               && isequal( varargin{2}, varargin{5} ) ... % compare Y-coordinates 
               && isequal( varargin{3}, varargin{6} ) ;   % compare Z-coordinates 
    end

    if ~exist( 'isSame' )
        error( 'See HELP MrdiGrid.comparegrids' ) ;
    end

end

