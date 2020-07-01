function [ basis ]= spherical_harmonics( orders, X, Y, Z )
%spherical_harmonics Return orthonormal spherical harmonic basis set 
% 
% SYNTAX
%
%    basis = spherical_harmonics( orders, X, Y, Z )
% 
% DESCRIPTION
%
% Returns a double-precision array of SH basis fields with the order/degree
% index along the 4th dimension.
%
% INPUTS
%
%    orders
%      Degrees of the desired terms in the series expansion, specified as a
%      vector of non-negative integers (`[0:1:n]` yields harmonics up to n-th
%      order)
%
%    X, Y, Z
%      Identically-sized 2- or 3-D numeric arrays of grid coordinates 
% 
% EXAMPLE
% 
% ```
% % Initialize grid positions
% [X,Y,Z] = ndgrid([-10:10],[-10:10],[-10:10]); 
%
% % 0th-to-2nd order terms inclusive
% orders = [0:2];
%
% basis = spherical_harmonics(orders, X, Y, Z);
% ```
%
% - `basis(:,:,:,1)` corresponds to the 0th-order constant term
%
% - `basis(:,:,:,2:4)` to 1st-order linear terms
%   - 2: *y*  
%   - 3: *z*   
%   - 4: *x*   
%
% - `basis(:,:,:,5:9)` to 2nd-order terms
%   - 5: *xy*  
%   - 6: *zy*   
%   - 7: *z2*   
%   - 8: *zx* 
%   - 9: *x2y2*  
%  
% NOTES 
%
% Based on calc_spherical_harmonics_arb_points_cz.m by jaystock@nmr.mgh.harvard.edu
    arguments
        orders(1,:) {mustBeNumeric,mustBeNonnegative};
        X {mustBeNumeric};
        Y {mustBeNumeric};
        Z {mustBeNumeric};
    end

%% Check inputs
assert( isequal(ndims(X), ndims(Y), ndims(Z)), ...
    'Input arrays X, Y, and Z must be identically sized') ;

switch ndims(X) 
    case 2
        gridSize = [ size(X) 1 ] ;
    case 3
        gridSize = size(X) ;
    otherwise 
        error('Input arrays X, Y, and Z must have 2 or 3 dimensions') ;
end

%% Initialize variables
nVoxels  = numel(X);
nOrders  = numel(orders) ;

harm_all = zeros( nVoxels, 1 ) ;

ii=0;

%% Compute basis 
for iOrder = 1 : nOrders 
   
    n = orders(iOrder);
    m = -orders(iOrder):1:orders(iOrder);

    for mm = 1 : numel(m)

        ii = ii+1;

        % The first 2*n(1)+1 columns of the output correspond to harmonics of
        % order n(1), and the next 2*n(2)+1 columns correspond to harmonics of
        % order n(2), etc. 
        harm_all(:,ii) = leg_rec_harmonic_cz( n, m(mm), X(:), Y(:), Z(:));
    end

end

%% Reshape to initial gridSize
nBasis = size( harm_all, 2) ;
basis  = zeros( [gridSize nBasis] ) ;

for iBasis = 1 : nBasis 
    basis(:,:,:, iBasis) = reshape( harm_all(:, iBasis), gridSize ) ;
end

% ----------------
%% Local functions

function out = leg_rec_harmonic_cz(n, m, pos_x, pos_y, pos_z)
% returns harmonic field for the required solid harmonic addressed by 
% n, m based on the iterative Legendre polynomial calculation
% Positive m values correspond to cosine component and negative to sine
%
% returned fields will eventually follow RRI's convention
% pos_... can be both value and vector/matrix

    r2        = pos_x.^2+pos_y.^2+pos_z.^2;
    r         = r2.^0.5;
    phi       = atan2(pos_y, pos_x);
    cos_theta = cos(atan2((pos_x.^2+pos_y.^2).^0.5, pos_z));
    %cos_theta=pos_z./r;

    if m>=0,
        c=1;
    else
        c=0;
        m=-m;
    end

    Ymn      = leg_rec(n,m,cos_theta);

    rri_norm = factorial(n+m+1)/factorial(n-m)/ffactorial(2*m);

    out      = (n+m+1)*r.^(n).*(cos(m*phi)*c+sin(m*phi)*(1-c)).*Ymn/rri_norm;

function out = ffactorial(n)
%FFACTORIAL FFactorial (double factorial) function.

    N = n(:);
    if any(fix(N) ~= N) || any(N < 0) || ~isa(N,'double') || ~isreal(N)
      error('MATLAB:factorial:NNegativeInt', ...
            'N must be a matrix of non-negative integers.')
    end

    if n==0 || n==1
        out=1;
    else
        out=n*ffactorial(n-2);
    end

end

function out=leg_rec(n, m, u)
% compute legendre polynomial values for dat using recursive relations

    if m>0
        p_mm=(-1)^m*ffactorial(2*m-1)*(1-u.^2).^(m/2);
    else
        p_mm=1;
    end

    if (n==m)
        out=p_mm;
    else
        p_mm1=(2*m+1)*u.*p_mm;
        
        if (n==m+1)
            out=p_mm1;
        else
            % recursive calculation needed
            a=m+2;
            p_ma_2=p_mm;
            p_ma_1=p_mm1;
            
            while 1
                p_ma=((2*a-1)*u.*p_ma_1-(a+m-1)*p_ma_2)/(a-m);
                
                if a==n
                    break;
                end
                % prepare next iteration
                p_ma_2=p_ma_1;
                p_ma_1=p_ma;
                a=a+1;
            end
            
            out=p_ma;
        end
    end
end %leg_rec

end %leg_rec_harmonic_cz

end % main
