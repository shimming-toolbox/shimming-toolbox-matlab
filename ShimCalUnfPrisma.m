% =========================================================================
% function [] = acquirecalibrationdata20170323( )

baselineDacValues.x    = 98.13 ;
baselineDacValues.y    = -1172.15 ;
baselineDacValues.z    = -1495.16 ; % note: dicom filenames were mislabeled as  _A01 instead of _A10
baselineDacValues.z2   = 4.61 ;
baselineDacValues.zx   = 0.47 ;
baselineDacValues.zy   = 5.15 ;
baselineDacValues.x2y2 = -63.91 ;
baselineDacValues.xy   = -21.38 ;



dacValues = [ baselineDacValues.x, (baselineDacValues.x - 30.0)   , (baselineDacValues.x + 30.0) ;
              baselineDacValues.y, (baselineDacValues.y - 30.0)   , (baselineDacValues.y + 30.0) ;
              baselineDacValues.z, (baselineDacValues.z - 30.0)   , (baselineDacValues.z + 30.0) ;
              baselineDacValues.z2, (baselineDacValues.z2 - 600.0) , (baselineDacValues.z2 + 600.0) ;
              baselineDacValues.zx, (baselineDacValues.zx - 600.0)  , (baselineDacValues.zx + 600.0) ;
              baselineDacValues.zy, (baselineDacValues.zy - 600.0)  , (baselineDacValues.zy + 600.0) ;
              baselineDacValues.x2y2, (baselineDacValues.x2y2 - 800.0)  , (baselineDacValues.x2y2 + 800.0) ;
              baselineDacValues.xy, (baselineDacValues.xy - 800.0)  , (baselineDacValues.xy + 800.0) ; ] ;
