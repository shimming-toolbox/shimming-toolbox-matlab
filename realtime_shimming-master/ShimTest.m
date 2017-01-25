classdef ShimTest < ShimUse
%SHIMTEST - Shim Testing
%
% .......
% 
% Description
%
%   ShimTest contains various methods for testing the 24-channel shim.
%
% =========================================================================
% Notes
% 
% ShimTest is a ShimUse subclass [ShimTest < ShimUse]
%
% Part of series of classes pertaining to shimming:
%
%     ProbeTracking
%     ShimCal
%     ShimCom
%     ShimOpt
%     ShimSpecs
%     ShimUse
%     ShimTest 
%
% =========================================================================
% Updated::20160827::ryan.topfer@polymtl.ca
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
%
% Usage
%
% elapsedTime = TESTRESPONSETIME( Shim, command )
% elapsedTime = TESTRESPONSETIME( Shim, command, nCycles )
% 
% Inputs
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
function [] = shimsine( Shim, Params )
%SHIMSINE 
%
% Issue sinusoidal current waveform   

Params.maxCurrent       = 0.25 ; % [units: A]
Params.nTestCycles      = 100 ;
Params.period           = 1 ; % period of waveform [units: s]
Params.updateFrequency  = 1/0.25 ; % update freq of shims [units: Hz]
Params.testTime         = 60 ;

% channelsToBankKey = Shim.getchanneltobankkey ;
% activeChannelIndices = channelsToBankKey(:,4) ;


% assert( (testChannel>0) && testChannel<=Shim.Specs.nActiveChannels) )


% if ~exist(Params.testTime 
Params.testTime = Params.nTestCycles * Params.period ;

currents = zeros(24,1);

figure(1)
axis([0 Params.testTime -0.51 0.51])
ylabel('Current (A)');
xlabel('Time (s)');

t=0;
for iChannel = 1 : 24
    hold on;
    plot( t, currents(iChannel),'*' ) ;
end
drawnow;

tic 
while t < Params.testTime
    t = toc ;
    currents( : ) = Params.maxCurrent * sin( t*Params.updateFrequency  ) ;
    % Shim.setandloadallshims( currents )
    % Shims.setallshims( zeros(Shims.Specs.Amp.nChannels, 1) ) ;
    disp(currents(1));
    
    for iChannel = 1 : 24
        
        figure(1)
        hold on;
        plot( t, currents(iChannel),'*' ) ;
    end
    drawnow;
    pause(0.1);
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
end
