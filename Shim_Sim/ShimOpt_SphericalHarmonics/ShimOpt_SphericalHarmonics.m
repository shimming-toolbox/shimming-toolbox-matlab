classdef ShimOpt_SphericalHarmonics < ShimOpt
%SHIMOPTSHARMONICS - Shim optimization using spherical harmonic basis
% 
% ShimOpt_SphericalHarmonics is a ShimOpt subclass. See ShimOpt documentation
% for general usage.
%
% .......
% 
% Usage
%
%       Shim = ShimOpt_SphericalHarmonics( Field )
%       Shim = ShimOpt_SphericalHarmonics( Field, Params )
% 
%   As for other ShimOpt objects, Field is a MaRdI-type object representing the
%   Field to be shimmed. The only difference is the (optional) input struct
%   Params, for which 2 (mutually exclusive) fields are configurable:
%   
%   Params.ordersToGenerate : default = [1:2]
%       
%       If set, this generates arbitary/ideal spherical harmonics.
%
%       .ordersToGenerate is a linear sequence of non-negative integers
%       specifying the orders of spherical harmonics to be generated and placed
%       along the 4th dimension of the array Shim.img 
%
%       e.g. default of [1:2], Shim.img will have 8 shim terms (three
%       1st-order, five 2nd-order harmonics) with values defined at the
%       voxel positions of the input Field object. 
%
%       For more info, See doc
%       ShimOpt_SphericalHarmonics.generatebasisfields 
%
%
%   Params.systemName : default = []
%
%       If set to 'IUGM_Prisma_fit' or 'HGM_Prisma', the returned Shim object 
%       possesses analogous terms to the respective Prisma scanner, with identical 
%       ordering along the 4th dimension of Shim.img, i.e. increments along this
%       dimension correspond to X, Y, Z, Z2, ZX, ZY, X2-Y2, and XY terms
%       respectively. The difference between the resulting ShimOpt_SphericalHarmonics
%       object and the corresponding ShimOpt_x_PrismaX object is that the shim
%       terms of the former are *ideally* generated as opposed to empirically
%       mapped.
%
% .......
%
% NOTE
%   The essential method GENERATEBASISFIELDS() is based on
%   calc_spherical_harmonics_arb_points_cz.m by jaystock@nmr.mgh.harvard.edu
%     
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

properties  
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt_SphericalHarmonics( varargin )
%SHIMOPTSHARMONICS - Shim Optimization with spherical harmonic basis set

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;

[ Field, Params ] = ShimOpt.parseinput( varargin ) ;

Params = ShimOpt_SphericalHarmonics.assigndefaultparameters( Params ) ;

if myisfieldfilled( Params, 'systemName' )
    switch Params.systemName
        case 'IUGM_Prisma_fit'
            Shim.System.Specs = ShimSpecs_IUGM_Prisma_fit() ;
        case 'HGM_Prisma'
            Shim.System.Specs = ShimSpecs_HGM_Prisma() ;
        otherwise
            error('Unimplemented system: Params.systemName must either be "IUGM_Prisma_fit" or "HGM_Prisma"') ;
    end
else  
    % Define ShimSpecs for the virtual shim system:
    Specs.systemName           = [ 'SphericalHarmonics_' ...
        num2str(min(Params.ordersToGenerate)) '-' num2str(max(Params.ordersToGenerate)) ] ;

    Specs.nChannels            = 0 ;
    for iOrder = 1 : length( Params.ordersToGenerate )
        Specs.nChannels = Specs.nChannels + 2*Params.ordersToGenerate(iOrder) + 1 ;
    end
    Specs.nActiveChannels      = Specs.nChannels ;

    Specs.maxCurrentPerChannel = Inf*ones( Specs.nChannels, 1 ) ;
    Specs.staticChannels       = true( Specs.nChannels, 1 ) ;
    Specs.dynamicChannels      = true( Specs.nChannels, 1 ) ;

    Specs.channelNames = cell( Specs.nChannels, 1 ) ;
    Specs.channelUnits = cell( Specs.nChannels, 1 ) ;

    for iCh = 1 : Specs.nChannels
        Specs.channelNames(iCh) = { [ 'SH_' num2str(iCh) ]} ;
        Specs.channelNames(iCh) = { [ '[AU]' ]} ;
    end

    Shim.System.Specs = ShimSpecs_Sim( Specs ) ;
end

%% -----
if ~isempty( Field ) 
    Shim.setoriginalfield( Field, Params.ordersToGenerate ) ;
end

end
% =========================================================================
function [] = setoriginalfield( Shim, Field, ordersToGenerate )
%SETORIGINALFIELD 
%
% [] = SETORIGINALFIELD( Shim, Field, ordersToGenerate )
%
% Sets Shim.Field
%
% Field is a FieldEval type object with .img in Hz

if nargin < 3
    error('Not enough input arguments.') ;
end

Shim.Model.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 )  ;

Shim.Field = Field.copy() ;

[X,Y,Z]  = Field.getvoxelpositions();
dR = Field.isocenter() ;
% Voxel positons are in the patient coordinate system,
% shift along Z to account for possible table displacement:
Z  = Z - dR(3) ; 

switch Shim.System.Specs.Id.systemName
    case {'IUGM_Prisma_fit', 'HGM_Prisma'} 
        Shim.img = Shim.generatebasisfields_prisma( X, Y, Z ) ;
    otherwise
        Shim.img = ShimOpt_SphericalHarmonics.generatebasisfields( ordersToGenerate, X, Y, Z ) ;
end


Shim.Hdr = Field.Hdr;

Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;

end
% =========================================================================
function [Corrections] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% Corrections = OPTIMIZESHIMCURRENTS( Shim, Params )

if nargin < 2 
    Params.dummy = [];
end

Corrections = optimizeshimcurrents@ShimOpt( Shim, Params ) ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================
methods(Access=private)
% =========================================================================
function [ basisFields ] = generatebasisfields_prisma( Shim, X, Y, Z )
%GENERATEBASISFIELDS_PRISMA
% 
% Wraps to ShimOpt_SphericalHarmonics.generatebasisfields(), reorders, and
% rescales the basis set to correspond to the Prisma

assert( nargin == 4 )
basisFields  = [] ;
basisFields0 = ShimOpt_SphericalHarmonics.generatebasisfields( [1:2], X, Y, Z ) ;

%% ------
% Returned basisFields are ordered like: 
%   Y, Z, X, XY, ZY, Z2, ZX, X2-Y2
% Reorder them in line with Siemens shims: 
%   X, Y, Z, Z2, ZX, ZY, X2-Y2, XY

Gx = basisFields0(:,:,:,3) ;
Gy = basisFields0(:,:,:,1) ;
Gz = basisFields0(:,:,:,2) ;

SZ2   = basisFields0(:,:,:,6) ;
SZX   = basisFields0(:,:,:,7) ;
SZY   = basisFields0(:,:,:,5) ;
SX2Y2 = basisFields0(:,:,:,8) ;
SXY   = basisFields0(:,:,:,4) ;

%% ------
% Rescale to the nominal values of the Prisma shims given in the 3d shim adjustments card:
% i.e. X,Y,Z should correspond to 1 micro-T/m, and 2nd order terms should be in micro-T/m^2

dbstop in ShimOpt_SphericalHarmonics at 204 ;
[minX, iMinX ] = min(X(:)) ;
[maxX, iMaxX ] = max(X(:)) ;

[minY, iMinY ] = min(Y(:)) ;
[maxY, iMaxY ] = max(Y(:)) ;

[minZ, iMinZ ] = min(Z(:)) ;
[maxZ, iMaxZ ] = max(Z(:)) ;

Gx = Gx*(1E6)*( Gx(iMaxX) - Gx(iMinX) )/(maxX - minX) ;
Gy = Gy*(1E6)*( Gy(iMaxY) - Gy(iMinY) )/(maxY - minY) ;
Gz = Gz*(1E6)*( Gz(iMaxZ) - Gz(iMinZ) )/(maxZ - minZ) ;

%% ------
Interpolant = scatteredInterpolant() ;

Interpolant.Method              = 'linear' ;
Interpolant.ExtrapolationMethod = 'linear' ;

% The following avoids the error from scatteredInterpolant when one
% attempts to form a 3d interpolant from a 2d input: 
isValidDim0 = [ numel(unique(X(:))) numel(unique(Y(:))) numel(unique(Z(:))) ] > 1 ;
r0          = [X(:) Y(:) Z(:)] ;

Interpolant.Points = r0(:, isValidDim0) ;

Interpolant.Values = SZ2(:) ;

basisFields = basisFields0 ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [ basisFields ]= generatebasisfields( orders, X, Y, Z )
%GENERATEBASISFIELDS
%  
% Generates orthonormal spherical harmonic (SH) basis set 
%
% .......
%  
% Usage
%
% [ basisFields ] = GENERATEBASISFIELDS( orders, X, Y, Z )
%
% Returns array of SH basis fields where the 4th dimension is the order/degree index
% e.g. if orders = [0:2],
%
%       basisFields(:,:,:,1) corresponds to the 0th order term,
%
%       basisFields(:,:,:,2:4) to 1st order terms
%         2 -> (y)  
%         3 -> (z)   
%         4 -> (x)   
%
%       basisFields(:,:,:,5:9) to 2nd orders
%         5 -> (xy)  
%         6 -> (zy)   
%         7 -> (z2)   
%         8 -> (zx) 
%         9 -> (x2y2)  
%
%       etc.
%
% Input
%
%   orders 
%       vector of non-negative integers orders to calculate (e.g. [0:3]).
%       Typically orders is specified as 0:1:N to obtain spherical harmonics
%       up to order N
%
%   X, Y, Z
%       3d arrays specifying voxel coordinates at which to calculate the harmonics 
%  
% .......
%
% Based on calc_spherical_harmonics_arb_points_cz.m by jaystock@nmr.mgh.harvard.edu


% set the orders_to_calculate as a vector of integers n with 0 < n < N. This 
% calculates accompanying spherical harmonic orders.  The first 2*n(1)+1 columns 
% of the output correspond to harmonics of order n(1), and the next
% 2*n(2)+1 columns correspond to harmonics of order n(2), etc.  %
% note that Hetherington refers to n as "degree" instead of "order"
%

assert( all( orders >=0 ) ) ;

gridSize = size(X) ;
nVoxels  = numel(X);

nOrders  = numel(orders) ;


harm_all = zeros( nVoxels, 1 ) ;

ii=0;

for iOrder = 1 : nOrders 
   
    n = orders(iOrder);
    m = -orders(iOrder):1:orders(iOrder);

    for mm = 1 : numel(m)

        ii = ii+1;

        harm_all(:,ii) = leg_rec_harmonic_cz( n, m(mm), X(:), Y(:), Z(:));

    end

end

nBasisFields = size( harm_all, 2) ;
basisFields  = zeros( [gridSize nBasisFields] ) ;

for iBasisField = 1 : nBasisFields 

    tmpField = reshape( harm_all(:, iBasisField), gridSize ) ; 
    
    % % normalize each term
    % basisFields(:,:,:, iBasisField) = tmpField/max(tmpField(:)) ;
    basisFields(:,:,:, iBasisField) = tmpField ;
end

%% change order of basis fields to correspond to Siemens Prisma


function out = leg_rec_harmonic_cz(n, m, pos_x, pos_y, pos_z)
% returns harmonic field for the required solid harmonic addressed by 
% n, m based on the iterative Legendre polynomial calculation
% Positive m values correspond to cosine component and negative to sine
%
% returned fields will eventually follow RRI's convention
% pos_... can be both value and vector/matrix

    r2=pos_x.^2+pos_y.^2+pos_z.^2;
    r=r2.^0.5;
    phi=atan2(pos_y, pos_x);
    cos_theta=cos(atan2((pos_x.^2+pos_y.^2).^0.5, pos_z));
    %cos_theta=pos_z./r;

    if m>=0,
        c=1;
    else
        c=0;
        m=-m;
    end

    Ymn=leg_rec(n,m,cos_theta);

    rri_norm=factorial(n+m+1)/factorial(n-m)/ffactorial(2*m);

    out=(n+m+1)*r.^(n).*(cos(m*phi)*c+sin(m*phi)*(1-c)).*Ymn/rri_norm;

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
end

end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static=true, Hidden=true)
% =========================================================================
function  [ Params ] = assigndefaultparameters( Params )
%ASSIGNDEFAULTPARAMETERS  
% 
% Params = ASSIGNDEFAULTPARAMETERS( Params )
% 
% DEFAULT_ORDERSTOGENERATE = [1:2] ;

DEFAULT_ORDERSTOGENERATE = [1:2] ;

if myisfieldfilled( Params, 'systemName' )
    assert( any(strcmp( Params.systemName, {'IUGM_Prisma_fit', 'HGM_Prisma'} )), ...
            'Unimplemented system: Params.systemName must either be "IUGM_Prisma_fit" or "HGM_Prisma"') ;
    
    if myisfieldfilled( Params, 'ordersToGenerate' )
        assert( all( Params.ordersToGenerate == [1:2] ), ...
            ['Incompatible Params.ordersToGenerate: The Prisma possesses 1st and 2nd order shims. ' ...
             'To simulate a higher order system, set Params.systemName = []'] ) ;
    end
end

if ~myisfieldfilled( Params, 'ordersToGenerate' ) 
   Params.ordersToGenerate = DEFAULT_ORDERSTOGENERATE ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end


