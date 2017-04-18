function[ roiOut ] = dialter( roiIn, R )
%DILATER dilates a 3D binary mask by 'R'
%
%   DILATER convolves binary input image ROI with an ellipsoid defined by
%   radius (or radii) R to return an expanded ROI.
%
%   Syntax
%
%   dilatedRoi = DILATER( roi, R)
%
%   Description
%
%   dilatedRoi = RILATER(ROI,R) returns ROI dilated by R. If R is a single
%   scalar, every dimension is dilated by R. If R is a 3-component vector
%   [Rx Ry Rz], each dimension is then dilated by its corresponding R value.
%
% see also SHAVER( ) 
%
% =========================================================================
% Updated::20170210::ryan.topfer@polymtl.ca
% =========================================================================


if any( R <= 0 )
	
    roiOut = roiIn ;
	disp('Shaver: Radius should be a positive integer. Returning input array') 

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

    if numel(R) == 1
        R = [R R R] ;
    end

    gridSize = size( roiIn ) ;
    
    sphericalKernel     = createellipsoid( 2*R + [1 1 1] , R) ;
    nPointsKernel       = sum( sphericalKernel(:) ) ;
	sphericalKernel     = sphericalKernel/nPointsKernel  ;

    % convolve with the kernel
    convSize = gridSize + size(sphericalKernel) - [1 1 1]; % *linear convolution size
    tmp = real(ifftn(fftn(roiIn,convSize).*fftn(sphericalKernel,convSize)));
   
    % retain regions of overlap btwn kernel & roi 
    tmp = tmp >= 1/(1 + nPointsKernel) ;
	
    roiOut = croparray( tmp, gridSize ) ;

	if inputIsFloat
	    roiOut = double( roiOut ) ;
	end

end

