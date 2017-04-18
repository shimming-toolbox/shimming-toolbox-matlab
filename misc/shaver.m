function[ roiOut ] = shaver( roiIn, R )
%SHAVER shaves a 3D binary mask by 'R'
%
%   SHAVER convolves binary input image ROI with an ellipsoid defined by
%   radius (or radii) R to return a contracted version of ROI
%
%   Syntax
%
%   ROIshaved = SHAVER(ROI,R)
%
%   Description
%
%   ROIshaved = SHAVER(ROI,R) returns ROI eroded by R. If R is a single
%   scalar, every dimension is eroded by R. If R is a 3-component vector
%   [Rx Ry Rz], each dimension is then eroded by its corresponding R value.
%   
%
%   see also DILATER( )
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
        R = [R R R];
    end
    
    gridSize = size( roiIn ) ;
	
    sphericalKernel     = createellipsoid( 2*R + [1 1 1] , R) ;
    nPointsKernel       = sum( sphericalKernel(:) ) ;
	sphericalKernel     = sphericalKernel/nPointsKernel  ;
    
    % convolve with the kernel
    convSize = gridSize + size(sphericalKernel) - [1 1 1]; % *linear convolution size
    tmp = real(ifftn(fftn(roiIn,convSize).*fftn(sphericalKernel,convSize)));
   
    % erode regions of non or partial overlap btwn kernel & roi 
    tmp = tmp > (1 - 0.99/nPointsKernel) ;
  
    % extract central region 
    roiOut = croparray( tmp, gridSize ) ;
    
	if inputIsFloat
	    roiOut = double( roiOut ) ;
	end

end
