%EVAL_ACDC_74
% 
% evaluate acdc_74 'real-time' field maps

%% load 'normal breathing' (Nb) data
MagNb   = MaRdI('/Volumes/mri/unf/acdc/acdc_74/acdc_74/04-realtime_fieldmap_normalBreath') ;
PhaseNb = MaRdI('/Volumes/mri/unf/acdc/acdc_74/acdc_74/05-realtime_fieldmap_normalBreath') ;

FieldsNb = FieldEval.mapfield( MagNb, PhaseNb ) ;
  
ProbeNb = ProbeTracking('/Volumes/mri/unf/acdc/acdc_74/PMU/acdc_74_01.resp');

%% load 'deep breathing' (Db) data
MagDb    = MaRdI('/Volumes/mri/unf/acdc/acdc_74/acdc_74/06-realtime_fieldmap_deepBreath') ;
PhaseDb  = MaRdI('/Volumes/mri/unf/acdc/acdc_74/acdc_74/07-realtime_fieldmap_deepBreath') ;

FieldsDb = FieldEval.mapfield( MagDb, PhaseDb ) ;

ProbeDb = ProbeTracking( '/Volumes/mri/unf/acdc/acdc_74/PMU/acdc_74_02.resp' ) ;

%% link respiratory PMU recordings
FieldsNb.associateaux( ProbeNb ) ;   
FieldsDb.associateaux( ProbeDb ) ;   

%% fit field to PMU recording
FieldNb = FieldEval.modelfield( FieldsNb );
FieldDb = FieldEval.modelfield( FieldsDb );

figure

subplot(221)
imagesc( FieldNb.img ) ;
axis equal
caxis([-200 200])
colorbar
title('Normal Breathing: Static B0') ;

subplot(223)
imagesc( FieldNb.Model.Riro.img ) ;
axis equal
caxis([-5 5])
colorbar
title('Normal Breathing: RMS RIRO') ;

subplot(222)
imagesc( FieldDb.img ) ;
axis equal
caxis([-200 200])
colorbar
title('Deep Breathing: Static B0') ;

subplot(224)
imagesc( FieldDb.Model.Riro.img ) ;
axis equal
caxis([-20 20])
colorbar
title('Deep Breathing: RMS RIRO') ;

%% compare pmu & field times series at arbitrary voxel in SC around ~C7
iRow = 20 ;
iCol = 30 ;
timeNb = FieldsNb.Aux.Data.t - FieldsNb.Aux.Data.t(1) ;
timeDb = FieldsDb.Aux.Data.t - FieldsDb.Aux.Data.t(1) ;

figure

subplot(221)
plot( timeNb, FieldsNb.Aux.Data.p, 'r' ) ;
hold on
plot( timeNb, FieldNb.Aux.Data.p*ones(size(timeNb)), '--k' ) ;
axis tight;
xlabel('Time (ms)')
ylabel('Resp. (A.U.)') ;
title('Normal Breathing') ;

subplot(223)
plot( timeNb, squeeze( FieldsNb.img(iRow, iCol,1,:) ), 'b' ) ;
hold on
plot( timeNb, FieldNb.img(iRow, iCol)*ones(size(timeNb)), '--k' ) ;
axis tight;
xlabel('Time (ms)')
ylabel('B0 (Hz)') ;

subplot(222)
plot( timeDb, FieldsDb.Aux.Data.p, 'r' ) ;
hold on
plot( timeNb, FieldDb.Aux.Data.p*ones(size(timeDb)), '--k' ) ;
axis tight;
xlabel('Time (ms)')
ylabel('Resp. (A.U.)') ;
title('Deep Breathing') ;

subplot(224)
plot( timeDb, squeeze( FieldsDb.img(iRow, iCol,1,:) ), 'b' ) ;
hold on
plot( timeNb, FieldDb.img(iRow, iCol)*ones(size(timeDb)), '--k' ) ;
axis tight;
xlabel('Time (ms)')
ylabel('B0 (Hz)') ;



