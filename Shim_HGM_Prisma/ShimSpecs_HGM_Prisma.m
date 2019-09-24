classdef ShimSpecs_HGM_Prisma < ShimSpecs
%SHIMSPECS_HGM_Prisma
% 
% Shim System Specifications for Siemens shims on Prisma at Montreal General Hospital
%    
%    ShimSpecs_HGM_Prisma is a ShimSpecs subclass
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimSpecs_HGM_Prisma(  )
%SHIMSPECS - Shim System Specifications 

Shim.Id.systemName   = 'HGM_Prisma' ;
Shim.Id.channelNames = { 'X (A11)' ; 'Y (B11)' ; 'Z (A10)' ; 'Z2 (A20)' ; 'ZX (A21)' ; 'ZY (B21)' ; 'X2-Y2 (A22)' ; 'XY (B22)' } ;
Shim.Id.channelUnits = { '[micro-T/m]' ; '[micro-T/m]'; '[micro-T/m]'; 
    '[micro-T/m^2]' ; '[micro-T/m^2]' ; '[micro-T/m^2]' ; '[micro-T/m^2]' ; '[micro-T/m^2]' ; } ;

Shim.Com = [] ;  
Shim.Dac = [] ;

% channels 1-3 refer to the gradient offsets, 4-8 to the 2nd order shims
Shim.Amp.nChannels       = 8 ;  
Shim.Amp.nActiveChannels = 8 ;

Shim.Amp.staticChannels  = true( Shim.Amp.nActiveChannels, 1 ) ;  
Shim.Amp.dynamicChannels = false( Shim.Amp.nActiveChannels, 1 ) ;  

% NB: One can use the Siemens commandline AdjValidate tool to get all the values below:
Shim.Amp.maxCurrentPerChannel = [ 2299.941, 2299.931, 2299.852, 2479.504, 1775.645, 1751.650, 1775.645, 1743.651 ]' ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [ shimValues ] = converttomultipole( shimValues )
%CONVERTTOMULTIPOLE
% 
% shimValues = CONVERTTOMULTIPOLE( shimValues )
%
% Shim values stored in MrProt (private Siemens DICOM.Hdr) are in units of 
% DAC counts for the gradient offsets and in units of mA for the 2nd order shims.
% CONVERTTOMULTIPOLE uses the information given by the Siemens commandline tool
%   AdjValidate -shim -info
% to convert a vector of shim settings in those units into the "multipole" values
% which are used in the Siemens GUI display (i.e. Shim3d)

nChannels = numel( shimValues ) ;

MAX_SHIMS.mp  = [ 2299.941; 2299.931; 2299.852; 2479.504; 1775.645; 1751.650; 1775.645; 1743.651 ] ;
MAX_SHIMS.dcm = [ 14370; 14323; 14078; 4999; 4999; 4999; 4999; 4999 ] ;

if nChannels == 3 
    % input shimValues are gradient offsets [units : DAC counts]
    % output shimValues units : micro-T/m]
    
    shimValues = shimValues .* ( MAX_SHIMS.mp(1:3) ./ MAX_SHIMS.dcm(1:3) ) ;

    if( any( abs( shimValues ) > MAX_SHIMS.mp(1:3) ) )
        warning('Multipole values exceed known system limits. Check inputs and limits specified in ShimSpecs') ;
    end

elseif nChannels == 5
    % input shimValues are for the 2nd order shims [units : mA]
    % output shimValues units : micro-T/m^2]

    shimValues = shimValues .* ( MAX_SHIMS.mp(4:8) ./ MAX_SHIMS.dcm(4:8) ) ;
    
    if( any( abs( shimValues ) > MAX_SHIMS.mp(4:8) ) )
        warning('Multipole values exceed known system limits. Check inputs and limits specified in ShimSpecs') ;
    end

elseif nChannels == 8

    shimValues = shimValues .* ( MAX_SHIMS.mp ./ MAX_SHIMS.dcm ) ;
    
    if( any( abs( shimValues ) > MAX_SHIMS.mp ) )
        warning('Multipole values exceed known system limits. Check inputs and limits specified in ShimSpecs') ;
    end

else
    error('Input shimValues vector should have 3, 5, or 8 elements, respectively corresponding to gradient terms, 2nd order shims, and both') ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
