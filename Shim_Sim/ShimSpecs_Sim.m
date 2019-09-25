classdef ShimSpecs_Sim < ShimSpecs
%SHIMSPECS_DES
% 
% Shim System Specifications for a simulated coil design, to be used in
% conjunction with ShimOpt_Sim
%
% First used in Lopez Rios et al. Integrated AC / DC coil and dipole Tx array
% for 7T MRI of the spinal cord. In: Proc. 27th Annu. Meet. ISMRM, Montreal,
% Canada, 2019. Abstr. 0220.
% 
% ShimSpecs_Sim is a ShimSpecs subclass. See ShimSpecs documentation for usage.
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
