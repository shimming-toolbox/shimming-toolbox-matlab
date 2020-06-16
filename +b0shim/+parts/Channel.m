classdef Channel 
%b0shim.parts.Channel Defines a shim coil element
%        
%    self = b0shim.parts.Channel()
%
% `Channel` provides a template description of a single coil element. The
% properties of a `Channel` object-array thereby describe an array of coils.
%
% See also  
% [ b0shim.parts.Contents ](./Contents.md)  
% [ b0shim.Config ](../Config.md)  

properties

    % Shim channel ID as a string-scalar.
    %
    % `name` defines the row names when tabulating optimization results 
    % in `b0shim.Opt.optimizeshimcurrents` (i.e. for print display).
    %
    % __EXAMPLE__ 
    % 
    % % Generically, this could simply be  
    % `name = "Ch. 1"`; 
    %
    % % Or, for a 2nd-order spherical harmonic term  
    % `name = "XY (B22)"`; 
    name(1,1) string  = "Ch" ;
    
    % Shim channel physical units as a string-scalar.
    %
    % `units` is used when tabulating optimization results in
    % `b0shim.Opt.optimizeshimcurrents` (i.e. for print display).
    %
    % __EXAMPLE__ 
    % 
    % % Units are often expressed in Amperes  
    % `units = "A"`;
    %
    % % Alternatively, for a 2nd-order spherical harmonic shim, units could be
    % `units = "mT/m^2"`;
    units(1,1) string = "A" ;
    
    % Range of permitted values as 2-element vector.
    % 
    % `limits` defines the default constraints during shim optimization:
    % `limits(1)` specifies the minimum value; 
    % `limits(2)` specifies the maximum value.
    % 
    % Both elements should be given in the corresponding `units`.
    %
    % __EXAMPLE__
    % 
    % For an "AC/DC" multi-coil array, with `units="A"`, the value might be 
    % `limits = [-2.5 2.5];`. 
    limits(1,2) {mustBeNumeric,mustBeReal} ; 

    % Label specifying coil positioning mode as a string-scalar.
    % 
    % Options for `positioning` are:  
    %
    % - "iso-fixed": The coil is immobile relative to isocenter (e.g. a gradient coil)
    %
    % - "table-fixed": The coil position changes with the patient-table (e.g. spine shim array)
    % 
    % - "": Otherwise (e.g. coils can be freely rearranged)
    positioning(1,1) string {mustBeMember(positioning,["iso-fixed" "table-fixed" ""])} = "iso-fixed";
end

% =========================================================================
% =========================================================================
methods 
% =========================================================================
function self = Channel( props )
    return
end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
