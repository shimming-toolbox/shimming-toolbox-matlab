function[ roiOut ] = dialter( roiIn, R )
%DILATER dilates a 3D binary mask by 'R'
%
%   DILATER convolves binary input image ROI with an ellipsoid defined by
%   radius (or radii) R to return an expanded ROI
% =========================================================================
%
%   Syntax
%
%   dilatedRoi = DILATER( roi, R)
%
%   Description
%
% Also see:
%   shaver.m
%
% =========================================================================
% Created: RT, topfer@ualberta.ca
%   Mon  9 Nov 2015 16:23:26 EST
% =========================================================================

if any( R == 0 )
	roiOut = roiIn ;
	disp('Radius cannot be zero. Returning input array') 
else
	inputIsFloat = true;
    tmp = whos('roiIn');
	switch tmp.class
	    case {'single','double'}
	        inputIsFloat   = true ;
	        roiIn          = logical( roiIn ) ;
	    case {'logical'}
	        inputIsFloat = false ;
	end

	if any( mod( size( roiIn ), 2 ) == 0 )
	    isCroppingToOddDimensions = true ;
	    roiIn                     = makeodd( roiIn ) ;
	    gridDimensionVector       = size( roiIn ) ;
	else
	    isCroppingToOddDimensions = false ;
	    gridDimensionVector = size( roiIn ) ;
	end

	sphere              = createellipsoid( gridDimensionVector, R) ;
	roiOut              = ifftc( fftc(roiIn) .* fftc( sphere/sum( sphere(:)) ) ) ;
	roiOut              = abs(roiOut) >= ( 1/( 1 + sum(sphere(:) ) ) ) ;

	if inputIsFloat
	    roiOut = double( roiOut ) ;
	end

	if isCroppingToOddDimensions
	    roiOut = makeodd( roiOut, 'isUndoing' ) ;
	end

end

