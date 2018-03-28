% ACDC SHIM DEMO
%
%
% =========================================================================
% Updated::20180328::ryan.topfer@polymtl.ca
% =========================================================================

cd('/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/scripts')

Shims = ShimUse( declareshimparamsdemo_201803() ) ;


% -------
% TODO 
% TEST CONNECTION 

% =========================================================================
% ACQUIRE + LOAD + PROCESS TRAINING DATA (PRESSURE TRACE + GRE FIELD MAPS)
% =========================================================================
Shims.acquiretrainingdata(): 

% -------
% TODO 
% ------
% delete this once above is tested
Shims.Data.Aux.Tracker = cell(1,1,3) ;
TmpTracked = ProbeTracked( ) ;
[TmpTracked.Data.p] = TmpTracked.loadmeasurementlog( [Shims.Params.dataLoadDir '20171013T150203-pressureLog-INS.bin'] ) ;
Shims.Data.Aux.Tracker{1} = TmpTracked ;

TmpTracked = ProbeTracked( ) ;
[TmpTracked.Data.p] = TmpTracked.loadmeasurementlog( [Shims.Params.dataLoadDir '20171013T150341-pressureLog-EXP.bin'] ) ;
Shims.Data.Aux.Tracker{2} = TmpTracked ;

TmpTracked = ProbeTracked( ) ;
[TmpTracked.Data.p] = TmpTracked.loadmeasurementlog( [Shims.Params.dataLoadDir '20171013T144940-pressureLog-BreathingTest.bin'] ) ;
Shims.Data.Aux.Tracker{3} = TmpTracked ;
% -------


Shims.loadandprocesstrainingdata() ; % process field maps (FieldEval-type) now in Shims.Data.Img{1,2,:}

% =========================================================================
% PREP OPTIMIZATION
% =========================================================================
% by default, shim VOI should begin as... Shims.Params.cordVoi ?

% delete ?
Params.shimVoi = Shims.Params.cordVoi ;
Params.shimVoi(1:32,:,:) = 0 ;

Shims.Opt.setshimvolumeofinterest( Params.shimVoi ) ;

% Shims.Opt.setshimvolumeofinterest( ...
%     Shims.Opt.getvaliditymask( { Shims.Opt.Field ; Shims.Opt.Field.Model.Shift } ) ) ;
%
% nii( Shims.Opt.Field.Hdr.MaskingImage ) ;
%
% % alternatively, avg. insp./exp echoes first? may be more robust...
% Params.dataWeights1 = Shims.Opt.derivedataweights( ...
%     { Shims.Data.Img{ 1, 1, 1 }  ; Shims.Data.Img{ 2, 1, 1 } }, 16, ...
%     Shims.Params.cordVoi ) ;
%
% Params.dataWeights2 = Shims.Opt.derivedataweights( ...
%     { Shims.Data.Img{ 1, 1, 2 }  ; Shims.Data.Img{ 2, 1, 2 } }, 16, ...
%     Shims.Params.cordVoi ) ;
%
% Params.dataWeights = ( Params.dataWeights1 + Params.dataWeights2 )/2 ;
%
% Params.dataWeights(1:32,:,:) = 0 ; % ignore brain...
%
% Shims.Params.dataWeights = Params.dataWeights ;
%
% Shims.Model.dataWeights

% =========================================================================
% OPTIMIZATION
% =========================================================================
% Shims.Params.dataWeights = Shims.Params.cordVoi ;
% Shims.Params.dataWeights(1:32,:,:) = 0 ;
% Shims.Opt.setshimvolumeofinterest( Params.dataWeights ) ;
Shims.Params.isReturningPseudoInverse = true ;

[currents] = Shims.Opt.optimizeshimcurrents( Shims.Params ) ;



PredictedFieldInspired =  Shims.predictshimmedfield( ) ;


predictedFieldInspired = Params.shimVoi .* PredictedFieldInspired.img ;
nii(predictedFieldInspired) ;



% % -------
% % Estimate optimal imaging frequencies for inspired condition:
% Shims.Opt.setoriginalfield( FieldInspired ) ;
% Shims.Opt.setshimvolumeofinterest( Params.shimVoi ) ;
% Shims.Model.currents = Params.Inspired.currents ;
%
% [f0, Params.Inspired.larmorShimsOff, Params.Inspired.larmorShimsOn] = ...
%     Shims.Opt.optimizelarmor( Params.shimVoi ) ;
%
%
% % -------
% % Estimate optimal imaging frequencies for expired condition:
% Shims.Opt.setoriginalfield( FieldExpired ) ;
% Shims.Opt.setshimvolumeofinterest( Params.shimVoi ) ;
%
% Shims.Model.currents = Params.Expired.currents ;
%
% [f0, Params.Expired.larmorShimsOff, Params.Expired.larmorShimsOn] = ...
%     Shims.Opt.optimizelarmor( Params.shimVoi ) ;
%
%
% % -------
% % Estimate optimal imaging frequencies for real-time acquisition:
%
% Params.larmorShimsOff = ...
%     ( Params.Inspired.larmorShimsOff + Params.Expired.larmorShimsOff )/2 ; 
%
% Params.larmorShimsOn = ...
%     ( Params.Inspired.larmorShimsOn + Params.Expired.larmorShimsOn )/2 ; 
%
% vpa(Params.larmorShimsOff)
%
% vpa(Params.larmorShimsOn)




% =========================================================================
% CALIBRATE REAL-TIME UPDATES 
% =========================================================================
Params = Shims.Opt.calibraterealtimeupdates( Params ) ;

% =========================================================================
% SET  
% =========================================================================

% TODO : PROPER CURRENT LIMITING IN REAL-TIME SHIMMING
Params.maxCurrents = max( [Params.Inspired.currents Params.Expired.currents], [],2) ;
Params.minCurrents = min( [Params.Inspired.currents Params.Expired.currents], [],2) ;

% Params.minClipThreshold        = Params.Expired.medianP ;
% Params.maxClipThreshold        = Params.Inspired.medianP ;
Params.minClipThreshold        = 9 ;
Params.maxClipThreshold        = 13 ;


% TODO : SAVE PARAMS and Shim.Opt 
save( [Params.dataLoadDir datestr(now,30) '-Params' ], 'Params' ) ;
% save( [Params.dataLoadDir datestr(now,30) '-ShimOpt' ], Shims.Opt ) ;

% =========================================================================
% SET SHIMS
% =========================================================================

% % -------

% -------
% Inspired : Real-time shimming (RTS) 
Params.measurementLogFilename         = [Params.dataLoadDir datestr(now,30) '-pressureLog-Inspired-RTS-GRE.bin'] ;
Params.sampleTimesFilename            = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Inspired-RTS-GRE.bin'] ;
Params.updateTimesFilename            = [Params.dataLoadDir datestr(now,30) '-updateTimes-Inspired-RTS-GRE.bin'] ;
Shims.runrealtimeshim( Params ) ;

% -------
% Inspired : Real-time shimming (RTS)
Params.measurementLogFilename         = [Params.dataLoadDir datestr(now,30) '-pressureLog-Inspired-RTS-EPI.bin'] ;
Params.sampleTimesFilename            = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Inspired-RTS-EPI.bin'] ;
Params.updateTimesFilename            = [Params.dataLoadDir datestr(now,30) '-updateTimes-Inspired-RTS-EPI.bin'] ;
Shims.runrealtimeshim( Params ) ;

% Expired : Real-time shimming (RTS)
Params.measurementLogFilename         = [Params.dataLoadDir datestr(now,30) '-pressureLog-Expired-RTS-GRE.bin'] ;
Params.sampleTimesFilename            = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Expired-RTS-GRE.bin'] ;
Params.updateTimesFilename            = [Params.dataLoadDir datestr(now,30) '-updateTimes-Expired-RTS-GRE.bin'] ;
Shims.runrealtimeshim( Params ) ;

% % -------
% Expired : Real-time shimming (RTS)
Params.measurementLogFilename         = [Params.dataLoadDir datestr(now,30) '-pressureLog-Expired-RTS-EPI.bin'] ;
Params.sampleTimesFilename            = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Expired-RTS-EPI.bin'] ;
Params.updateTimesFilename            = [Params.dataLoadDir datestr(now,30) '-updateTimes-Expired-RTS-EPI.bin'] ;
Shims.runrealtimeshim( Params ) ;


% -------
% Normal breathing RTS
Params.measurementLogFilename     = [Params.dataLoadDir datestr(now,30) '-pressureLog-Breathing-RTS.bin'] ;
Params.sampleTimesFilename        = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Breathing-RTS.bin'] ;
Params.updateTimesFilename        = [Params.dataLoadDir datestr(now,30) '-updateTimes-Breathing-RTS.bin'] ;
Shims.runrealtimeshim( Params ) ;

Params.measurementLogFilename     = [Params.dataLoadDir datestr(now,30) '-pressureLog-Breathing-RTS.bin'] ;
Params.sampleTimesFilename        = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Breathing-RTS.bin'] ;
Params.updateTimesFilename        = [Params.dataLoadDir datestr(now,30) '-updateTimes-Breathing-RTS.bin'] ;
Shims.runrealtimeshim( Params ) ;

% -------
% Normal breathing spine shim OFF
Params.pressureLogFilename = [Params.dataLoadDir datestr(now,30) '-pressureLog-Breathing-ShimOff.bin'] ;
Params.sampleTimesFilename = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Breathing-ShimOff.bin'] ;
[pressureLog, sampleTimes] =Shims.Opt.Tracker.recordandplotpressurelog( Params ) ;


% last inspired field map: static currents.
Shims.Com.setandloadallshims( Params.Inspired.currents ) 

% =========================================================================
% TURN OFF SHIMS
% =========================================================================
Shims.Com.resetallshims()


















