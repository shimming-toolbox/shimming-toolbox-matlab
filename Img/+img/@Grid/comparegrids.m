function isSame = comparegrids( varargin )
%COMPAREGRIDS Return `true` if voxel positions coincide 
%        
%     isSame = comparegrids( Grid1, Grid2 )
%     isSame = comparegrids( xyz1, xyz2 )
%     isSame = comparegrids( X1, Y1, Z1, X2, Y2, Z2 )
% 
% Wrapper function to `img.Grid.isequal()`
%
% __ETC__
%
% See also 
% img.Grid.isequal

    if nargin == 2
        if isa( varargin{1}, 'img.Grid' ) && isa( varargin{2}, 'img.Grid' )
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
        error( 'See HELP Grid.comparegrids' ) ;
    end

end

