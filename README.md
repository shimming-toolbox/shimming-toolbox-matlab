# shimming_rri

=========================================================================
SHIMMING_RRI
=========================================================================

This library consists of programs to perform shimming (static and real-time)
using the Resonance Research Inc. 24-channel spine shim [REF]

Series of classes pertaining to shimming:

    ProbeTracking
    ShimCal
    ShimCom
    ShimOpt
    ShimSpecs
    ShimUse
    ShimTest 

Type 'doc [class name]' for the documentation.


-------------------------------------------------------------------------
Usage:
    
    classdef ShimUse
        Shimming is controlled by the user through an object of type ShimUse.
    ShimUse is a subclass of ShimCom and thus inherits its methods.

    classdef ShimCom
        ShimCom is responsible for communicating with the RRI equipment.

    Methods

-------------------------------------------------------------------------
Notes:

-------------------------------------------------------------------------
Updated::20160826::ryan.topfer@polymtl.ca


-------------------------------------------------------------------------
Todo:

