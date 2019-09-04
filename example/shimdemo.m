% Example script to perform offline shim optimization based on Siemens gradient
% echo (time-series) data

%% -----
% load gradient echo time series:
Mag    = MaRdI( '~/Matlab/shimming/example/acdc_78p/GRE_FM_RTSHIM_0129' ) ;
Phase  = MaRdI( '~/Matlab/shimming/example/acdc_78p/GRE_FM_RTSHIM_0130' ) ;

% compute field maps (time-series):
Fields = FieldEval.mapfield( Mag, Phase ) ;
% note: there are unwrapping errors present in the field map at the top of the bottle due to the large voxels

%% -----
% load respiratory data
Probe = ProbeTracking( '~/Matlab/shimming/example/acdc_78p/20190825T202139_pressureProbeRec.mat' ) ;

% NOTE: Starting point for adapting code to make field-fitting + shimming comptatible with the Siemens PMU:
%
%   Loading this file will add struct called Data to the Matlab workspace
%
%   load('~/Matlab/shimming/example/acdc_78p/20190825T202139_pressureProbeRec.mat')
% 
%   Required is a script that takes the Siemens PMU recording as input and returns a file/struct equivalent to this. 
%
%
%   Essential are the Data fields:
%   
%   Data.pRaw    : raw signal time-series (a vector of doubles)
%   
%   Data.p       : processed signal (guessing this is all Siemens saves - in which case, Data.p and Data.pRaw can simply be copies)
%
%   Data.t       : time points associated with each signal value     [units: ms]
%
%   Data.trigger : vector whose entries are equal to 1 when Data.t featured a trigger, and 0 otherwise
%                  more precisely, Data.trigger(i) = 1 when a trigger occurred *within* the preceding time interval Data.t( i-1 : i )

%% -----
% link respiratory recording to field maps:
%
% Fields.associateaux( Probe, Params ) returns Field.Aux : 
% a respiratory recording object (like Probe itself) but with its corresponding
% Data. subfields now corresponding to the image times
%
% (i.e. Probe.Data is cropped and interpolated to form Field.Aux.Data )

% estimated transmission delay: 
% NOTE: for Siemens PMU, Params.auxDelay would presumably be = 0 [the default in ProbeTracking.associateaux()]
Params.auxDelay = 0.05 ; % [units: s] 

Fields.associateaux( Probe, Params ) ;   

%% -----
% Perform fitting of B0 to Aux resp signal:
%
% FieldEval.modelfield( Fields ) will return the fitted/modelled Field:
%
% Field.img :
%   the static B0 estimate
%
% Field.Model.Riro.img : 
%   the respiration-varying component
%   note, RIRO values are scaled by the root-mean square of the debiased respiratory signal 
%   i.e. rms( Field.Aux.Data.p - mean( Field.Aux.Data.p ) ) 

Field = FieldEval.modelfield( Fields ) ;


%% -----
% Create a ShimOpt object corresponding to the hardware
Params = [] ;
Shims  = ShimOpt_Greg( Params, Field ) ;

% The Shims.Aux contains the ShimOpt object of the Siemens Gradient+Shim coils (called: ShimOpt_IUGM_Prisma_fit), which is used for joint optimization.
% Note: If you instanciate Shims directly using using ShimOpt_IUGM_Prisma_fit, then the .Aux field will be empty.
% 
% Define the target region for shimming with a binary mask (same size as the Matlab variable: Field.img)
%
% by default, it will be the intersection of Field.Hdr.MaskingImage (where
% reliable field values exist) and Shims.getshimsupport() (where shim reference
% map values exist)
%
% if the user assigns a different region, only the overlap with the reliable
% region will be considered in the optimization.
shimVoi = Shims.getshimsupport() & Shims.Aux.getshimsupport() & Shims.Field.Hdr.MaskingImage ;
Shims.setshimvolumeofinterest( shimVoi ) ;

% Define which terms are to be included in the static shim optimization (boolean vector):
%   [TX_FREQ ; CUSTOM_COIL ; SIEMENS_COIL ]
% Where:
%   - TX_FREQ is the Tx frequency of the MR sytem (obtained via frequency adjustement routine)
%   - CUSTOM_COIL: Vector comprising each shim element. The name of each element is stored in: Shims.System.Specs.Id.channelNames  
%   - SIEMENS_COIL: Vector comprising the gradient+shim elements. The name of each element is stored in: Shims.Aux.System.Specs.Id.channelNames
% 
% NOTE: 
% If the "Shims" object is instantiated directly using the Siemens coil (Shims = ShimOpt_IUGM_Prisma_fit())
% then the .Aux does not exist, and the vector becomes: 
%   [TX_FREQ ; SIEMENS_COIL ]

% Activate real-time shimming
Params.isRealtimeShimming = 1 ;
% Define which shim terms are to be included in the real-time shim:
% e.g.
%   Tx-freq adjustment ; none of the 8 multi-coil channels ; Prisma gradients, but not the five 2nd order shims: 
Params.activeDynamicChannelsMask = [ true ; false(8, 1) ; true(3, 1); false(5, 1) ] ;

%% -----
% finally, perform the optimization:
% 
% Shims.optimizeshimcurrents() returns struct Currents with fields 
%
% Currents.static
%   static shim settings 
%   (Prisma values are in multi-pole units: micro-T/m for 1st order, micro-T/m^2 for 2nd order)
%
% Currents.realtime
%   'coupling coefficients' for real-time shimming 
%   (i.e. for a given respiratory measurement p, the corresponding shim currents for the real-time
%   correction would be p*Currents.realtime
% 
% vector entries are ordered the same way as the ChannelMask terms above
%
% static and realtime values are also respectively stored in
% Shims.Model.currents and Shims.Model.couplingCoefficients 
% as well as writen to .txt files
%
% Params.mediaSaveDir [default='./'] assigns the write directory for the results files 
Currents = Shims.optimizeshimcurrents( Params ) ; 

%NOTE: 
% The data here is not a great example insofar as all the shim reference maps
% were acquired with the "Greg" coil while the demo gre images were acquired
% with the head coil, which sits much closer to the table surface.
% (i.e. the reference maps only cover top part of the image volume)
%
% TODO: reacquire Prisma reference maps using large water container at ISO,
% using spine coil + flex coils on top
%
% furthermore there are unwrapping errors present in the field map at the top of the bottle
