This folder contains template classes ("abstract" classes in MATLAB
nomenclature) for a generalized shim system. Idea being that corresponding
subclasses will be defined for each specific shim system, but these
subclasses should all inherit from the following parent/template classes:

**ShimSpecs( )**

*System specifications Re: amplifier, DAC, etc.*

**ShimCom( )**

*Low-level communication with the amplifier/hardware.*

**ShimOpt( )**

*Shim optimization. Paths to shim reference maps (.mat) files need to be specified here, and system-specific features (e.g. nonlinear current constraints) may need to be accounted.*

For more information, see the documentation for each class, e.g. in the MATLAB
prompt, type HELP ShimOpt()
