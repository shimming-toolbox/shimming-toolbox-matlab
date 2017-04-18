classdef ShimOptSHarmonics < ShimOpt
%SHIMOPTSHARMONICS - Shim optimization using spherical harmonic basis
%
% .......
% 
% Usage
%
% Shim = ShimOptSHarmonics( Params, Field )
% 
% Field is a MaRdI-type object representing the Field to be shimmed.
%
% Defaults 
%
% Params.isGeneratingBasis       = true ; 
% 
% Params.ordersToGenerate        = [0:2];
% 
% Params.pathToShimReferenceMaps = [] ; % if provided, then
%   Params.isGeneratingBasis = false (spherical harmonic basis set is not
%   generated)
% 
% Params.isInterpolatingReferencemaps = true ; % if true, the reference maps
%   are automatically interpolated to the grid (voxel positions) of 'Field'
%
% Params.TrackerSpecs = [] ;
%   .TrackerSpecs is a parameters struct for ProbeTracking(). See HELP
%   ProbeTracking() for more information.
%
%   Shim contains fields
%
%       .img
%           Shim reference maps
%
%       .Hdr
%           Info Re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%       .Field
%           Object of type MaRdI pertaining to field distribution to be shimmed
%
%       .Model
%           .currents  
%               Optimal shim current vector (i)
%               [units A]
%           .field     
%               Optimal shim field from projection of i onto reference maps (Ai)
%               [units Hz]
%           .couplingCoefficients
%               For realtime shimming, relates vector relating field to pressure
%               [units Hz/Pa]
%           .dcCurrentsOffsets
%               For realtime shimming, vector of "y-intercept" currents 
%               (i.e. currents for pressure = 0 Pa)
%               [units A]
%
%       .Tracker
%           Object of type Tracking (e.g. ProbeTracking() ) 
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
%
%    Tracking
%    ShimCal
%    ShimCom
%    ShimEval
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%
% ShimOpt is a MaRdI subclass [ShimOpt < MaRdI]
%     
% =========================================================================
% Updated::20170412::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
% .....
% OPTIMIZESHIMCURRENTS()
%   Run (unconstrained) Cg solver;
%   If solution achievable given system constraints, return solution;
%   Else, run fmincon given constraints & return that solution instead;
%
% =========================================================================

% properties % defined in parent class ShimOpt 
    % Field ; % object of type MaRdI
    % Model ;
    % Tracker ; % object of type ProbeTracking
% end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOptSHarmonics( Params, Field )
%SHIMOPTSHARMONICS - Shim Optimization with spherical harmonic basis set
%
% Params.isGeneratingBasis
% Params.ordersToGenerate

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params = ShimOptSHarmonics.assigndefaultparameters( Params ) ;

Shim.Model = [ ] ; 
Shim.Tracker = ProbeTracking( Params.TrackerSpecs )  ; 

% % .......
% % Load shim basis if provided 
% if ~isempty(Params.pathToShimReferenceMaps)
%     
%     Params.isGeneratingBasisSet = false 
%
%     ShimUse.display(['\n Preparing for shim ...  \n\n'...
%             'Loading shim reference maps from ' Params.pathToShimReferenceMaps '\n\n']) ;
%
%     load( Params.pathToShimReferenceMaps ) ;
%
% end


if Params.isGeneratingBasis || Params.isInterpolatingReferenceMaps

    assert( (nargin > 1) && ~isempty(Field), 'Must input Field [MaRdi-type object]. See HELP')
    
    Shim = Shim.setoriginalfield( Field ) ;

    if Params.isGeneratingBasis 

            [X,Y,Z]  = Field.getvoxelpositions();
            Shim.img = ShimOptSHarmonics.generatebasisfields( Params.ordersToGenerate, X, Y, Z );
            Shim.Hdr = Field.Hdr;

    elseif Params.isInterpolatingReferenceMaps
             
            Shim = interpolatetoimggrid( Shim, Field )
            Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;

    end

end


end
% =========================================================================
function Shim = optimizeshimcurrents( Shim, CgParams )
%OPTIMIZESHIMCURRENTS 
%
% Shim = OPTIMIZESHIMCURRENTS( Shim )
% Shim = OPTIMIZESHIMCURRENTS( Shim, CgParams )
%
% CgParams : Parameters struct for conjugate-gradient optimization
%
%   .tolerance     [default = 1E-6] 
%   .maxIterations [default = 10000]

if nargin < 1
    error('Function requires at least 1 argument of type ShimOpt')
elseif nargin == 1
    CgParams.dummy = [];
end

% Params for conjugate-gradient optimization
CgParams.tolerance     = 1E-6 ;
CgParams.maxIterations = 10000 ;    

A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

b = M*(-Shim.Field.img(:)) ;

% ------- 
% Least-squares solution via conjugate gradients
Shim.Model.currents = cgls( A'*M'*M*A, ... % least squares operator
                            A'*M'*b, ... % effective solution vector
                            zeros( [Shim.getnactivechannels() 1] ), ... % initial model (currents) guess
                            CgParams ) ;
    
Shim = Shim.setforwardmodelfield ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================
methods(Access=protected)
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
%   i.e. 
%       basisFields(:,:,:,1) corresponds to 0th order, 
%
%       basisFields(:,:,:,2:4) to 1st orders
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

    nTermsithOrder = numel(m);

    for mm = 1 : nTermsithOrder
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
% Add default parameters fields to Params without replacing values (unless empty)
%
% DEFAULT_PATHTOSHIMREFERENCEMAPS = [] ;
% DEFAULT_PROBESPECS = [] ;
%
% DEFAULT_ISINTERPOLATINGREFERENCEMAPS = true ;
%
% DEFAULT_ISGENERATINGBASIS       = true ;
% DEFAULT_ORDERSTOGENERATE        = [0:2] ;

DEFAULT_PATHTOSHIMREFERENCEMAPS = [] ;
DEFAULT_TRACKERSPECS            = [] ;

DEFAULT_ISINTERPOLATINGREFERENCEMAPS = true ;

DEFAULT_ISGENERATINGBASIS       = true ;
DEFAULT_ORDERSTOGENERATE        = [0:2] ;

if ~myisfield( Params, 'pathToShimReferenceMaps' ) || isempty(Params.pathToShimReferenceMaps)
   Params.pathToShimReferenceMaps = DEFAULT_PATHTOSHIMREFERENCEMAPS ;
end

if ~myisfield( Params, 'TrackerSpecs' ) || isempty(Params.TrackerSpecs)
   Params.TrackerSpecs = DEFAULT_TRACKERSPECS ;
end

if ~myisfield( Params, 'isInterpolatingReferenceMaps' ) || isempty(Params.isInterpolatingReferenceMaps)
   Params.isInterpolatingReferenceMaps = DEFAULT_ISINTERPOLATINGREFERENCEMAPS ;
end

if ~myisfield( Params, 'isGeneratingBasis' ) || isempty(Params.isGeneratingBasis)
   Params.isGeneratingBasis = DEFAULT_ISGENERATINGBASIS ;
end

if ~myisfield( Params, 'ordersToGenerate' ) || isempty(Params.ordersToGenerate)
   Params.ordersToGenerate = DEFAULT_ORDERSTOGENERATE ;
end

end
% =========================================================================

end

end


