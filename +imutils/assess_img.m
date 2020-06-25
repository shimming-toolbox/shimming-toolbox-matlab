function stats = assess_img( img, gridSpacing, voi, filename )
%ASSESS_IMG Return
%
%     stats = assess_img( img )
%     stats = assess_img( img, gridSpacing, voi )
%     stats = assess_img( img, gridSpacing, voi, filename )
    arguments
        img(:,:,:,1) double {mustBeNumeric};
        gridSpacing(1,3) double {mustBeNumeric};
        voi(:,:,:,1) {valid.mustBeBoolean} = isfinite(size(img));
        filename {valid.mustBeStringOrChar} = tempname;
    end
% 
% voi 
%    binary array the same size as Field.img indicating the region of
%    interest over which field calculations are made. 
%    default: Field.Hdr.MaskingImage
%
% filename
%   output to text file using writetable()
%
% stats contains fields
%
%   .volume
%       volume of region of interest (voi) [units: cm^3]
%
%   .mean
%       mean value of the field over the voi
%
%   .median
%       median value of the field over the voi
%
%   .std
%       standard deviation of Field.img over the voi
%   
%   .rmsePerCm3
%       L2 norm of the field (i.e. residual) over the voi normalized by the volume
%
%   .meanAbs
%       mean absolute value of the field over the voi.
%
%   .medianAbs
%       median absolute value of the field over the voi
%
%   .min
%   .max
assert( ndims( img ) <=3, 'Multiple volumes not supported. TODO' )

if nargin < 2 || isempty(voi)
    voi = isfinite( size(img) ) ;

end

voi = logical( voi ) ;

stats.volume     = nnz( voi ) .* prod( 0.1*Field.getvoxelspacing() )  ; % [units: cm^3]
stats.mean       = mean( img( voi ) ) ;
stats.median     = median( img( voi ) ) ;
stats.std        = std( img( voi ) ) ;
stats.rmsePerCm3 = norm( img( voi ), 2 )/stats.volume ;
stats.meanAbs    = mean( abs( img( voi ) ) ) ;
stats.stdAbs     = std( abs( img( voi ) ) ) ;
stats.min        = min( ( img( voi ) ) ) ;
stats.max        = max( ( img( voi ) ) ) ;

if nargin == 3 && ischar( filename ) 
    measure = {'Volume (cm^3)'; 'Mean (Hz)' ; 'Median (Hz)' ; 'St. dev. (Hz)' ; 'RMSE/Volume (Hz/cm^3)' ; 'Mean[abs.] (Hz)'; 'Std[abs.] (Hz)'; 'Min (Hz)'; 'Max (Hz)'} ;
    value   = num2str([ stats.volume ; stats.mean ; stats.median ; stats.std ; stats.rmsePerCm3 ; stats.meanAbs ; stats.stdAbs ; stats.min ; stats.max ;], 4 ) ;
    writetable( table( measure, value ), filename ) ;
end

end

