classdef ShimTest < ShimUse
%SHIMTEST
%
% .......
% 
%   ShimTest is a ShimUse subclass containing miscellaneous functions for
%   testing the 24-channel shim array.
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%     ProbeTracking
%     ShimCom
%     ShimOpt
%     ShimSpecs
%     ShimUse
%     ShimTest 
%
% =========================================================================
% Updated::20160802::ryan.topfer@polymtl.ca
% =========================================================================
properties   

end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimTest(  )
%SHIMTEST  

end
% =========================================================================
function [elapsedTime] = testresponsetime( Shim, command, nCycles )
%TESTRESPONSETIME
%
% Repeatedly queries MXD to test rapidity of command & response. 
% .......
%   Usage
%
% elapsedTime = TESTRESPONSETIME( Shim, command )
% elapsedTime = TESTRESPONSETIME( Shim, command, nCycles )
% 
% .......
%   Inputs
%
%   command
%       function handle (i.e. to a shim command)
%   
%   nCycles
%       number of test cycles
%
% .......
%
% e.g. Issue 'getsystemheartbeat' 100 times.
% 
% >> TESTRESPONSETIME( Shim, @Shim.getsystemheartbeat, 100 ) ;
%
% Prints:
% ??????


DEFAULT_NCYCLES   = 1000 ;

if nargin < 3
    nCycles = DEFAULT_NCYCLES ;
end

tic;

for iCycle = 1 : nCycles ;
   command( ) ; 
end

elapsedTime = toc ;
Shim.display(['Elapsed time: ' num2str(elapsedTime) ' s.' ]) ;
Shim.display(['N commands issued: ' num2str(nCycles) ]) ;
Shim.display(['Pause between commands: ' num2str(Params.pauseTime) ' s.' ]) ;
Shim.display(['Avg. time per command (accounting for pauses)\n: ' ...
    num2str( (elapsedTime-nCycles*Params.pauseTime)/nCycles) ' s.' ]) ;

end
% =========================================================================
function [] = shimsine( Shim, testChannel, Params )
%SHIMSINE
%   

% Params.currentAmplitude = 1 ; % [units: A]
% Params.nTestCycles      = 10 ;
% Params.period           = 1 ; % period of waveform [units: s]
% Params.updateFrequency  = 1/0.05 ; % update freq of shims [units: Hz]
%
% assert( (testChannel>0) && testChannel<=Shim.Specs.nActiveChannels) )
%
% channelsToBankKey = Shim.getchanneltobankkey ;
%
% testTIme = Params.nTestCycles * Params.period ;
%
% t=0;
% tic 
% % systemResponse = Shim.setandloadshim( bankIndex, channelIndexByBank, current ) 
%     
% while t < testTime
%     t = toc ;
%     currents( testChannel ) = currentAmplitude sin( t*Params.updateFrequency  ) ;
%
%     % % write updateallshims() function!
%     % for iChannel = 1 : Shim.Specs.nActiveChannels
%         ShimsC.setandloadshim( channelsToBankKey(iChannel, 2), channelsToBankKey(iChannel,3), currents( iChannel ) ) ;
%     % end
% end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
end
