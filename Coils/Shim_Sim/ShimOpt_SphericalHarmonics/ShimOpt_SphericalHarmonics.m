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
        Shim.img = ShimOpt_SphericalHarmonics.generatebasisfields_siemens( X, Y, Z ) ;
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
methods(Static)
% =========================================================================
function [ basisFields ]= generatebasisfields( orders, X, Y, Z )

basisFields = b0shim.compute.spherical_harmonics( orders, X, Y, Z);

end
% =========================================================================
function [ basisFields ] = generatebasisfields_siemens( X, Y, Z )

basisFields = b0shim.compute.siemens_basis([1:2],X,Y,Z);

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


