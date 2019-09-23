classdef ShimOpt_SphericalHarmonics < ShimOpt
%SHIMOPTSHARMONICS - Shim optimization using spherical harmonic basis
%
% .......
% 
% Usage
%
% Shim = ShimOpt_SphericalHarmonics( Params, Field )
% 
% Field is a MaRdI-type object representing the Field to be shimmed.
%
% Defaults 
% 
% Params.ordersToGenerate        = [1:2];
%
%   Shim contains fields
% ... TODO : doc
% =========================================================================
% 
% ShimOpt_SphericalHarmonics is a ShimOpt subclass. See ShimOpt documentation
% for general usage.
%     
% =========================================================================
% Author::ryan.topfer@polymtl.ca
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

properties  
    Params;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt_SphericalHarmonics( Params, Field )
%SHIMOPTSHARMONICS - Shim Optimization with spherical harmonic basis set
% 
% Params.ordersToGenerate

Shim.img   = [] ;
Shim.Hdr   = [] ;
Shim.Field = [] ;       
Shim.Model = [] ;
Shim.Aux   = [] ;

Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

[ Field, Params ] = ShimOpt.parseinput( varargin ) ;

Shim.Params = ShimOpt_SphericalHarmonics.assigndefaultparameters( Params ) ;

% =========================================================================
% Define ShimSpecs for virtual shim system:
% =========================================================================
Specs.systemName           = [ 'SphericalHarmonics_' ...
    num2str(min(Shim.Params.ordersToGenerate)) '-' num2str(max(Shim.Params.ordersToGenerate)) ] ;

Specs.nChannels            = 0 ;
for iOrder = 1 : length( Shim.Params.ordersToGenerate )
    Specs.nChannels = Specs.nChannels + 2*Shim.Params.ordersToGenerate(iOrder) + 1 ;
end
Specs.nActiveChannels      = Specs.nChannels ;

Specs.maxCurrentPerChannel = Inf*ones( Specs.nChannels, 1 ) ;
Specs.staticChannels       = true( Specs.nChannels, 1 ) ;
Specs.dynamicChannels      = true( Specs.nChannels, 1 ) ;

Specs.channelNames = cell( Specs.nChannels, 1 ) ;

for iCh = 1 : Specs.nChannels
    Specs.channelNames(iCh) = { [ 'SH_' num2str(iCh) ]} ;
end

Shim.System.Specs = ShimSpecs_Des( Specs ) ;



if ~isempty( Field ) 
    Shim.setoriginalfield( Field ) ;
end

end
% =========================================================================
function [] = setoriginalfield( Shim, Field, currents )
%SETORIGINALFIELD 
%
% [] = SETORIGINALFIELD( Shim, Field )
% [] = SETORIGINALFIELD( Shim, Field, currents )
%
% Sets Shim.Field
%
% Field is a FieldEval type object with .img in Hz

if nargin < 2
    error('Not enough input arguments.') ;
elseif nargin == 2
    currents = 0;
    warning('Assuming field map was acquired with all shim channels at 0 A.');
end

Shim.Model.currents = currents ;
Shim.Field = Field.copy() ;


[X,Y,Z]  = Field.getvoxelpositions();
dR = Field.getisocenter() ;
% Voxel positons are in the patient coordinate system,
% shift along Z to account for possible table displacement:
Z  = Z - dR(3) ; 
Shim.img = ShimOpt_SphericalHarmonics.generatebasisfields( Shim.Params.ordersToGenerate, X, Y, Z );
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
% function currents = optimizeshimcurrents( Shim, CgParams )
% %OPTIMIZESHIMCURRENTS 
% %
% % Shim = OPTIMIZESHIMCURRENTS( Shim )
% % Shim = OPTIMIZESHIMCURRENTS( Shim, CgParams )
% %
% % CgParams : Parameters struct for conjugate-gradient optimization
% %
% %   .tolerance     [default = 1E-6] 
% %   .maxIterations [default = 10000]
%
% if nargin < 1
%     error('Function requires at least 1 argument of type ShimOpt')
% elseif nargin == 1
%     CgParams.dummy = [];
% end
%
% % Params for conjugate-gradient optimization
% CgParams.tolerance     = 1E-10 ;
% CgParams.maxIterations = 100000 ;    
%
% A = Shim.getshimoperator ;
% M = Shim.gettruncationoperator ;
%
% b = M*(-Shim.Field.img(:)) ;
%
% % ------- 
% % Least-squares solution via conjugate gradients
% Shim.Model.currents = cgls( A'*M'*M*A, ... % least squares operator
%                             A'*M'*b, ... % effective solution vector
%                             zeros( [Shim.getnactivechannels() 1] ), ... % initial model (currents) guess
%                             CgParams ) ;
%     
% currents = Shim.Model.currents ;
%
% end
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

DEFAULT_ORDERSTOGENERATE        = [1:2] ;

if ~myisfield( Params, 'ordersToGenerate' ) || isempty(Params.ordersToGenerate)
   Params.ordersToGenerate = DEFAULT_ORDERSTOGENERATE ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end


