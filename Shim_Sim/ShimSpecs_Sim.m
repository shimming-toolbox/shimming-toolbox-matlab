classdef ShimSpecs_Sim < ShimSpecs
%SHIMSPECS_SIM
% 
% Shim system specifications for a simulated coil design. 
% 
% Usage
%
%   Specs = ShimSpecs_Sim( Specs )
%
% **This class only serves to typecast a parameters struct defined in a ShimOpt_
% class (e.g. ShimOpt_SphericalHarmonics, or ShimOpt_Sim) into a "ShimSpecs"
% object given that certain functions may expect the latter as opposed to a
% generic struct.**
%
% First used in Lopez Rios et al. Integrated AC / DC coil and dipole Tx array
% for 7T MRI of the spinal cord. In: Proc. 27th Annu. Meet. ISMRM, Montreal,
% Canada, 2019. Abstr. 0220.
% 
% ShimSpecs_Sim is a ShimSpecs subclass. For more info, see ShimSpecs.m
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

properties
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimSpecs_Sim( Specs )
%SHIMSPECS - Shim System Specifications 

Shim.Com = [] ; % N/A
Shim.Dac = [] ; % N/A

Shim.Id.systemName            = Specs.systemName ;
Shim.Id.channelNames          = Specs.channelNames ;
Shim.Id.channelUnits          = Specs.channelUnits ;

Shim.Amp.nChannels            = Specs.nChannels ;
Shim.Amp.nActiveChannels      = Specs.nChannels ;
Shim.Amp.maxCurrentPerChannel = Specs.maxCurrentPerChannel ;
Shim.Amp.staticChannels       = Specs.staticChannels ;
Shim.Amp.dynamicChannels      = Specs.dynamicChannels ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
