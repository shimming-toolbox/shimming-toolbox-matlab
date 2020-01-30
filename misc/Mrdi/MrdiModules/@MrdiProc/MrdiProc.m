classdef (Abstract) MrdiProc < handle
%MrdiProc   MR DICOM Image Processing 
%
% Member methods for image processing 
% 
% e.g.
%   filter()
%   resliceimg()
%   setmaskingimage()
%   ...etc.

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Proc = MrdiProc(  )

end
% =========================================================================
% =========================================================================
end

% =========================================================================
% =========================================================================    
methods
    %.....
    [] = filter( Img, weights, Params )
    %.....
    [F] = regrid( Img, X_Ep, Y_Ep, Z_Ep, varargin )
    %.....
    [] = setmaskingimage( Img, mask )
    %.....
end
% =========================================================================
% =========================================================================
methods(Access=protected)
    %.....
    Interpolant = getinterpolant( Img, method, extrapolationMethod )
end
% =========================================================================
% =========================================================================
methods(Hidden=true)
% NOTE
%   Temporarily placing cropimg() and zeropad() here since these methods
%   1) might be deprecated
%   2) may or may not be useful
    %.....
    Img = cropimg( Img, gridSizeImgCropped, centralPoint )
    %.....
    Img = zeropad( Img, padSize, padDirection )
end
% =========================================================================
% =========================================================================
methods
% =========================================================================
function timeAverage = timeaverage( Img )
%TIMEAVERAGE  Return mean( Img.img, 5 ) 
%
% Img = TIMEAVERAGE( Img ) 

    timeAverage = mean( Img.img, 5 ) ;

end
% =========================================================================
function timeStd = timestd( Img )
%TIMESTD Return std( Img.img, 0, 5 ) 
%
% standardDeviation = TIMESTD( Img ) 

    timeStd = std( Img.img, 0, 5 ) ;

end
% =========================================================================

% =========================================================================
% =========================================================================
end


end
