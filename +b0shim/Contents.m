%b0shim Shim system representation, configuration, and implementation
%
% ## Intro
%
% The **Shimming-Toolbox** has an ambitious aim: To be compatible with shim
% systems *generally*. This is a tall order as even *existing* devices vary
% considerably in their respective designs and implementations. While a
% *plug-'n-play* "out-of-the-(tool)box" solution is, frankly, out of the
% question, the hope is to produce software that works for the greatest range
% of systems with the least amount of wheel-reinvention (system-specific code)
% possible.
%
% To this end, the code-base employs templates and abstract classes, which is
% largely reflected in the specific directory structure of the `b0shim`
% package. Hence, a "big picture" overview of the package structure/logic is
% outlined here and geared toward developers. Individual package components
% and "getting started" examples are described elsewhere. 
%
% ## Components
%
% ### Overview
%
% `b0shim` contains code and configuration files serving to describe and enable
% the use of b0-shim coils:
%
% ```
% +b0shim:
% |
% ├── @Opt: Shim optimization class 
% |
% ├── @Com: (Abstract) shim communication interface class
% |
% ├── @Config: Shim system configuration template class
% |
% └── +parts: Shim system-component template classes 
% |       |
% |       ├── @Channel: A shim coil-element
% |       |
% |       └── @Port: Serial port parameters (shim microcontroller communication) 
% |
% └── +coils: System-specific subpackages (json configuration files and custom subclasses)
%         |  
%         ├── +greg: 8-channel AC/DC neck coil
%         |       |
%         |       ├── config.json
%         |       |
%         |       └── @Com
%         |
%         ├── +rriyan: 24-channel spine array
%         |       |
%         |       ├── config.json
%         |       |
%         |       ├── @Opt
%         |       |
%         |       └── @Com  
%         |
%         ├── +IUGM_Prisma_fit: Siemens 3T system at UNF
%         |       |
%         |       ├── config.json
%         |       |
%         |       └── @Opt
%         |
%         └── + [Etc.]
% ```
%
% At the top-level of the package, there are three classes folders:
%
% 1. `@Opt`: A generic class for handling the actual shim optimization. 
% **NOTE**: WIP, to be refactored!
%
% 2. `@Com`: An abstract class defining the USB-serial communication interface
% for custom multi-channel coils. (The class does not apply to MRI scanner
% shims or virtual systems.)
%
% 3. `@Config`: Defines the template and serves as the standard datatype for
% shim system configurations. 
%
% Also, two subpackage folders:
%
% 1. `+parts`: Contains component classes (e.g. `@Channel, @Port`), which are
% used to define the general system configuration, i.e. the properties of
% `@Config` are *in part* composed of these objects).
%
% 2. `+coils`: Contains the system-specific subpackages, i.e. the code here is
% the actual *implementation* of a particular shim system.
%
% ### Integration 
% 
% The code-base represents specific systems using two basic components, each
% with a standard datatype (MATLAB class) and a corresponding file convention
% for saving to disk:
%
%| Component                                          |   `class`     |*file* | 
%|:---------------------------------------------------|:-------------:|:-----:|
%|Configuration parameters (e.g. max. current/channel)|`b0shim.Config`|*.json*| 
%|Reference maps (images of Hz-shift/unit-current)    |`b0shim.Basis` |*.nii* |
%
% __TODO__ The Ref. image component is a work in progress: currently formed
% by the `ShimOpt` class and stored as `.mat`, but this is slated to change!


% A `b0shim.Config` object is a required argument for both `b0shim.Opt` and `b0shim.Com`.
%
%
% namely, to *do* anything with a given shim system—optimize it for a field
% map, or control the hardware via MATLAB (for applicable systems), the
% configuration parameters need to be defined.
%
%
% For descriptions of the parameters involved in the system configuration, see
% the documentation for [`b0shim.Config`](./Config.md).
% 
% ... To be continued!


% ```mermaid
% graph TD 
%     start[("b0shim")]
%     start --> opt("@Opt")
%     start --> com("@Com")
%     start --> specs("@Specs")
%     start --> parts[("+parts")]
%     parts --> ch("@Channel")
%     parts --> port("@ComPort")
%     start --> coils[("+coils")];
%     coils --> coil1[("+greg")];
%     coils --> coil2[("+rriyan")];
%     coils --> coil3[("+IUGM_Prisma_fit")];
%     coil1 --> c11["specs.json"];
%     coil1 --> c12("@Com");
%     coil2 --> c21["specs.json"];
%     coil2 --> c22("@Com");
%     coil2 --> c23("@Opt");
%     coil3 --> c31("specs.json");
%     coil3 --> c32("@Opt");
%     specs -.-|"composes^"|ch;
%     specs -.-|"composes^"|port;
%  ```
%
%
