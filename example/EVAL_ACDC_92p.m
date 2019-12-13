%EVAL_ACDC_92p
% 

pathToMag   = '/Users/evaalonsoortiz/Documents/Academic/Postdoc_Julien/Projects/Dynamic_Shimming/Data/acdc_92p/06-realtime_fieldmap_1slice_TR6p8_4x4x4mm' ;
pathToPhase = '/Users/evaalonsoortiz/Documents/Academic/Postdoc_Julien/Projects/Dynamic_Shimming/Data/acdc_92p/07-realtime_fieldmap_1slice_TR6p8_4x4x4mm' ;

Fields = FieldEval( pathToMag, pathToPhase ) ; 

Pmu   = ProbeTracking( '/Users/evaalonsoortiz/Documents/Academic/Postdoc_Julien/Projects/Dynamic_Shimming/Data/acdc_92p/06-realtime_fieldmap_1slice_TR6p8_4x4x4mm/acdc_92p_realtime_fieldmap_1slice_TR6p8_4x4x4mm.resp' ) ;


%% link respiratory PMU recordings
Fields.associateaux( Pmu ) ;  

%% fit field to PMU recording
Field = FieldEval.modelfield( Fields ) ;

figure

subplot(2,1,1);
imagesc( Field.img ) ;
axis equal
caxis([-200 200])
colorbar
title('Static B0') ;

subplot(2,1,2);
imagesc( Field.Model.Riro.img ) ;
axis equal
caxis([-1 1])
colorbar
title('RMS RIRO') ;



%% compare pmu & field times series at arbitrary voxel
iRow = 95 ;
iCol = 37 ;
time = Fields.Aux.Data.t - Fields.Aux.Data.t(1) ;


figure

subplot(2,1,1);
plot( time, Fields.Aux.Data.p, 'r' ) ;
hold on
plot( time, Field.Aux.Data.p*ones(size(time)), '--k' ) ;
axis tight;
xlabel('Time (ms)')
ylabel('Resp. (A.U.)') ;
title('PMU') ;

subplot(2,1,2);
plot( time, squeeze( Fields.img(iRow, iCol,1,:) ), 'b' ) ;
hold on
plot( time, Field.img(iRow, iCol)*ones(size(time)), '--k' ) ;
axis tight;
xlabel('Time (ms)')
ylabel('B0 (Hz)') ;





