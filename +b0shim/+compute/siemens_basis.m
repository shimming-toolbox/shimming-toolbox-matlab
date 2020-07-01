function [ basis ] = siemens_basis( ~, X, Y, Z )
%siemens_basis Return 1st & 2nd order *ideal* Siemens shim fields 
%
% SYNTAX 
%    
%    basis = siemens_basis( orders, X, Y, Z )
% 
% DESCRIPTION
% 
% The function first wraps `b0shim.compute.spherical_harmonics()` to generate
% 1st and 2nd order spherical harmonic `basis` fields at the grid positions
% given by arrays `X,Y,Z`. *Following Siemens convention*, `basis` is then:
%
% - Reordered along the 4th dimension as *X, Y, Z, Z2, ZX, ZY, X2-Y2, XY*
%
% - Rescaled to Hz/unit-shim, where "unit-shim" refers to the measure displayed
% in the Adjustments card of the Syngo console UI, namely:
%
%   - 1 micro-T/m for *X,Y,Z* gradients (= 0.042576 Hz/mm)
%   - 1 micro-T/m^2 for 2nd order terms (= 0.000042576 Hz/mm^2)
%
% The returned `basis` is thereby in the form of ideal "shim reference maps",
% ready for optimization.
%
% INPUTS
%
%    orders (uint row-vector)
%      Degrees of the desired terms in the series expansion, specified as a
%      vector of non-negative integers (`[0:1:n]` yields harmonics up to n-th
%      order)
%
%    X (numeric 2- or 3-d array)
%      "Right->Left" grid coordinates in the patient coordinate system
%      (i.e. DICOM reference, units of mm)
%
%    Y (numeric 2- or 3-d array)
%      "Anterior->Posterior" grid coordinates in the patient coordinate system
%      (i.e. DICOM reference, units of mm)
%
%    Z (numeric 2- or 3-d array)
%      "Inferior->Superior" grid coordinates in the patient coordinate system
%      (i.e. DICOM reference, units of mm)
% 
% OUTPUTS
%
%    basis (double 4-D array)
%      Spherical harmonic basis fields
%
% NOTES/TODO
%
% 1. For now, `orders` is, in fact, ignored: fixed as [1:2]—which is suitable
% for the Prisma (presumably other Siemens systems as well)—however, the
% 3rd-order shims of the Terra should ultimately be accommodated too. (Requires
% checking the Adjustments/Shim card to see what the corresponding terms and
% values actually are). So, for now, `basis` will always be returned with 8 terms
% along the 4th dim.
%
% See also
% [ b0shim.compute.spherical_harmonics ]
    arguments
        ~ ;%orders(1,:) {mustBeNumeric,mustBeNonnegative};
        X {mustBeNumeric};
        Y {mustBeNumeric};
        Z {mustBeNumeric};
    end

sh = b0shim.compute.spherical_harmonics( [1:2], X, Y, Z ) ;

%% Reorder terms along 4th array dim. in line with Siemens shims: 
% X, Y, Z, Z2, ZX, ZY, X2-Y2, XY
sh = reordertosiemens( sh ) ; 

scalingFactors = computenormalizationfactors() ;
basis          = zeros( size( sh ) ) ; 

for iCh = 1 : size( sh, 4 ) 
   basis(:,:,:,iCh) = scalingFactors(iCh) * sh(:,:,:,iCh) ; 
end

return;

% ----------------
%% Local functions

function [ sh1 ] = reordertosiemens( sh0 )
%REORDERTOSIEMENS 
% 
%   sh1 = reordertosiemens( sh0 )
%
% Reorder 1st-2nd order basis terms along 4th dim. from
%
% 1. Y, Z, X, XY, ZY, Z2, ZX, X2-Y2 (output by b0shim.compute.spherical_harmonics), to
%
% 2. X, Y, Z, Z2, ZX, ZY, X2-Y2, XY (in line with Siemens shims)

assert( ( nargin == 1 ) && ( size( sh0, 4 ) == 8 ) )

sh1(:,:,:,1) = sh0(:,:,:,3) ;
sh1(:,:,:,2) = sh0(:,:,:,1) ;
sh1(:,:,:,3) = sh0(:,:,:,2) ;
sh1(:,:,:,4) = sh0(:,:,:,6) ;
sh1(:,:,:,5) = sh0(:,:,:,7) ;
sh1(:,:,:,6) = sh0(:,:,:,5) ;
sh1(:,:,:,7) = sh0(:,:,:,8) ;
sh1(:,:,:,8) = sh0(:,:,:,4) ;

end %reordertosiemens()

function [ scalingFactors ] = computenormalizationfactors()
%COMPUTENORMALIZATIONFACTORS
%
%  scalingFactors = computenormalizationfactors()
%
%  returns a vector of scalingFactors to apply to the (properly reordered)
%  ideal 1st+2nd order spherical harmonic fields to scale the terms as "shim
%  reference maps" in units of Hz/unit-shim 
% -----
% Gx, Gy, and Gz should yield 1 micro-T of field shift per metre
% equivalently, 0.042576 Hz/mm
%
% 2nd order terms should yield 1 micro-T of field shift per metre-squared
% equivalently, 0.000042576 Hz/mm^2
% 
% NOTE: The method has been worked out empirically 
% (further testing could be a good idea)

%% create basis on small 3x3x3 mm^3 isotropic grid
[XIso, YIso, ZIso] = meshgrid( [-1:1], [-1:1], [-1:1] ) ;

sh = b0shim.compute.spherical_harmonics( [1:2], XIso, YIso, ZIso ) ;

% Reorder terms along 4th array dim. in line with Siemens shims: 
sh = reordertosiemens( sh ) ; 

nChannels      = size( sh, 4) ; % = 8
scalingFactors = zeros( nChannels, 1 ) ;

% indices of reference positions for normalization: 
iX1   = find( ( XIso == 1 ) & ( YIso == 0 ) & ( ZIso == 0 ) ) ;
iY1   = find( ( XIso == 0 ) & ( YIso == 1 ) & ( ZIso == 0 ) ) ;
iZ1   = find( ( XIso == 0 ) & ( YIso == 0 ) & ( ZIso == 1 ) ) ;

iX1Z1 = find( ( XIso == 1 ) & ( YIso == 0 ) & ( ZIso == 1 ) ) ;
iY1Z1 = find( ( XIso == 0 ) & ( YIso == 1 ) & ( ZIso == 1 ) ) ;
iX1Y1 = find( ( XIso == 1 ) & ( YIso == 1 ) & ( ZIso == 0 ) ) ;

% order the reference indices like the sh field terms 
iRef = [iX1 iY1 iZ1 iZ1 iX1Z1 iY1Z1 iX1 iX1Y1]' ;

%% ------
% scaling:
% 1st order terms yield 1 micro-T of field shift per m (i.e 0.042576 Hz/mm )
% 2nd order terms yield 1 micro-T of field shift per m^2 (i.e 0.000042576 Hz/mm^2 )

% distance from iso/origin to adopted reference point [units: mm]
r = [1 1 1 1 sqrt(2) sqrt(2) 1 sqrt(2)] ;

% invert polarity of certain terms:
sh(:,:,:,[2,3,5,8]) = -sh(:,:,:,[2,3,5,8] ) ;

orders = [1 1 1 2 2 2 2 2] ;

for iCh = 1 : nChannels
    field = sh(:,:,:,iCh) ;
    scalingFactors(iCh) = 42.576*( ( r(iCh) * 0.001 )^orders(iCh) )/field( iRef( iCh ) ) ;
end

end %computenormalizationfactors()

end
