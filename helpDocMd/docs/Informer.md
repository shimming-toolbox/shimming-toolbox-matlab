# Informer

**Filetype:** _MATLAB&reg; classdef_

**Synopsis:** _Details the functionality of a .m file_

The `Informer` class informs the `Documentor` class about the contents of .m
files: providing the functional details needed to document them.

In general, the class works behind the scenes and would not be directly
called upon by a user; however, an Informer object can be constructed
independently as indicated below.

__CONSTRUCTOR SYNTAX__
      
      Info = Informer( mFile ) ;

Creates an `Informer` object pertaining to .m file (a script, function, class
method, or classdef file) pointed to by the path string `mFile`.
If `mFile` contains multiple files, then `Info` is returned as an
object-array.

`Informer` has but two properties: `mFile` and `Attributes`.

`Attributes` is a struct containing all available functional details
regarding `mFile`.

`Attributes` cannot be set directly, but is updated whenever `mFile` is set.

The read-only fields of `Attributes` depend on the given .m-file type and
should be fairly self-explanatory given the field names. More detail is
available in the method documentation for Informer.getmattributes.

    Documentation for Informer
       doc Informer


__ATTRIBUTES__


<table>
<table border=1><tr><th>Hidden</th><th>Sealed</th><th>Abstract</th><th>Enumeration</th><th>ConstructOnLoad</th><th>HandleCompatible</th><th>RestrictsSubclassing</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- InferiorClasses : [N/A] 
- ContainingPackage : [N/A] 
- EventList : [N/A] 
- EnumerationMemberList : [N/A] 
- SuperclassList : [N/A] 

- - -
## Properties


### mFile

**Synopsis:** _.m file path: Attributes property will update whenever mFile is set_

  .m file path: Attributes property will update whenever mFile is set

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : 
- SetMethod : Informer.set.mFile
- DefaultValue : /Users/ryan/Projects/General/scripts/shim/helpDocMd/src/@Informer/Informer.m
- Validation: 
Validator functions: mustBeFile
- DefiningClass : Informer

### Attributes

**Synopsis:** _Functional description of the .m file_

Attributes is a struct of .m file attributes, returned from a call to
Informer.getmattributes( mFile ). It contains the following basic
fields:

- mType: Type of .m file: string scalar returned from Informer.mfiletype( mFile ).
Possibilities are: ["script","function","classdef","method","NA"]
    
- .Name: Name of the script, function, class or class method

If the .m file is a valid MATLAB file (i.e. mType ~= "NA"), then
Attributes also contains fields:

- .Description: Header line of help-text (string vector returned from
    Informer.extracthelpheader) 

- .DetailedDescription: Body of help-text (string vector returned from
    Informer.extracthelpbody) 

### References ###

Remaining fields of Attributes vary depending on the type of file
(i.e. Attributes.mType). For more info, see the method documentation:
Informer.getmattributes

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- GetAccess : public
- SetAccess : private
- PartialMatchPriority : [N/A] 
- GetMethod : 
- SetMethod : 
- DefaultValue : 
- Validation: 
Class: struct
Validator functions: 
- DefiningClass : Informer

---
## Methods


---


### Informer

**Synopsis**: _Details the functionality of a .m file_ 

The `Informer` class informs the `Documentor` class about the contents of .m
files: providing the functional details needed to document them.

In general, the class works behind the scenes and would not be directly
called upon by a user; however, an Informer object can be constructed
independently as indicated below.

__CONSTRUCTOR SYNTAX__
      
      Info = Informer( mFile ) ;

Creates an `Informer` object pertaining to .m file (a script, function, class
method, or classdef file) pointed to by the path string `mFile`.
If `mFile` contains multiple files, then `Info` is returned as an
object-array.

`Informer` has but two properties: `mFile` and `Attributes`.

`Attributes` is a struct containing all available functional details
regarding `mFile`.

`Attributes` cannot be set directly, but is updated whenever `mFile` is set.

The read-only fields of `Attributes` depend on the given .m-file type and
should be fairly self-explanatory given the field names. More detail is
available in the method documentation for Informer.getmattributes.


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : mFile
- OutputNames : Info
- DefiningClass : Informer

---


### mfiletype

**Synopsis**: _Returns the type of .m file ["script"|"function"|"classdef"|"method"]_ 

     
     [mType, mPath, mExist] = mfiletype( mFile )

__DESCRIPTION__
`mfiletype` checks the list of Matlab source files specified by the [string-,
char-, or cellstr-] array of file paths `mFile` and returns:

1. `mType` a string vector with entries

`mType(i)`| when input `mFile(i)` is...
----------|-----------------------------------------------------------------|
"script"  | a script file
"function"| a function file (even one with void arguments)
"classdef"| a class definition file
"method"  | a .m function (~=constructor) in a folder beginning with "@"
    "NA"    | an unimplemented or non-Matlab file with a .m file extension
     ""     | an invalid file path (e.g. folder, non-Matlab, or non-existent)

2.`mPath` a string vector of the full file system paths.

3. `mExist` a vector of return values (doubles) from Matlab function `exist` [2],
namely: `mExist(i) = exist( mFile(i) );` the value *should* be == 2.

__NOTE__
Although `mType(i)` will be "NA" whenever file `i` could not be assessed
(e.g. a function that produces an error when called by `nargin`), an
assignment other than "NA" does *not* guarantee the file is in
working condition, and it may still contain errors.

This owes to `mfiletype`'s simplistic implementation [1]: The approach
relies entirely on MATLAB's `exist` function [2] which, as the official
documentation states, does not check the contents of '.m' files.

[1]: (https://blogs.mathworks.com/loren/2013/08/26/what-kind-of-matlab-file-is-this/)
[2]: (https://www.mathworks.com/help/matlab/ref/exist.html/)

See also
EXIST


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : mType, mPath, mExist
- DefiningClass : Informer

---


### getmattributes

**Synopsis**: _Return functional description of a .m file_ 

     
     Att = GETMATTRIBUTES( mFile )

Returns a struct of attributes `Att` describing the .m file `mFile` (a
script, function, classdef, or class method) given as a path string

__OUTPUT__

- Att: a struct of .m file attributes containing the following basic fields:

    - .mType: the type of .m file

    - .Name: Name of the script, function, class or class method 
        
    - .Description: Header line of help-text (string vector returned from
    Informer.extracthelpheader) 

    - .DetailedDescription: Body of help-text (string vector returned from
    Informer.extracthelpbody) 

1. Function files

If the .m file is a function, Att additionally contains fields:

- .nInputs: Number of input arguments

- .nOutputs: Number of output arguments

*TODO*: implement means of adding (at least for inputs): argument names,
types, validation functions, and defaults, when declared in an
arguments block, and add these as fields to Attributes)

2. Classdef files

If the .m file is a class definition, Attributes derives the following
additional fields from an instance of the associated meta.class object:
(For more info, refer to the MATLAB documentation: <https://www.mathworks.com/help/matlab/ref/meta.class.html>)

_Fields containing logical scalars_:

- .Hidden
- .Sealed
- .Abstract
- .Enumeration
- .ConstructOnLoad
- .HandleCompatible
- .RestrictsSubclassing

_Fields containing string arrays_:

- .SuperclassList (the names of superclasses)
- .InferiorClasses (the names of deriving inferior classes)
- .ContainingPackage (the name of the containing package, if applicable, as a string-scalar)

_Fields containing struct arrays_:

- .MethodList (class methods derived from meta.method objects)

    Elements of MethodList possess the following fields, with all but the final 4
    (which contain string arrays) containing scalar logicals:
    - .Static 
    - .Abstract 
    - .Static 
    - .ExplicitConversion
    - .Sealed
    - .Hidden
    - .Abstract
    - .Access 
    - .InputNames
    - .OutputNames
    - .DefiningClass 

- .PropertyList (Class properties derived from meta.property objects)

TODO: elaborate...


TODO: add events and enumerations substructs from meta.class:
- .EventList
- .EnumerationMemberList

#### Notes

Despite _sounding_ very useful based on the meta.class property names
(including those normally hidden but made visible upon converting a meta.class
object to a struct!), unfortunately, it seems that MathWorks hasn't yet
bothered to implement anything for many of the meta.class properties that
would actually make them meaningful (at least, not as of MATLAB 2019b).

(Hence, for now, although the .Description and .DetailedDescription fields of
the returned attributes struct borrow their field names from (nominal)
meta.class properties, these are fields are actually assigned independently
of meta.class by calls to Informer methods: gethelptext(), extracthelpheader()
and extracthelpbody().)

Furthermore, the current implementation of meta.class exhibits some curious behaviour:
e.g. it seems that if you add an arguments block to a class method---even if
you do not specify default values, which would normally render the arguments
optional---suddenly the InputNames of the corresponding meta.method object
disappear and are replaced by an uninformative 'varargin'??

__ETC__

See also

- Informer.mfiletype
- Informer.gethelptest
- Informer.gethelpbody
- Informer.gethelpheader
- <https://www.mathworks.com/help/matlab/ref/meta.class.html meta.class>
- <https://www.mathworks.com/help/matlab/ref/meta.property.html>
- <https://www.mathworks.com/help/matlab/ref/meta.validation-class.html meta.Validation>


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : Att
- DefiningClass : Informer

---


### gethelptext

**Synopsis**: _Return help-text of script|function|method|property|class as string-vector_ 


mHelp = GETHELPTEXT( name )


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : mHelp
- DefiningClass : Informer

---


### extracthelpheader

**Synopsis**: _Return the leading line of help-text_ 

     
      mHelpHeader = EXTRACTHELPHEADER( mHelp )
      mHelpHeader = EXTRACTHELPHEADER( mHelp, name )

If 'name' is provided as a second input argument and it appears as the first
word of the header line (irrespective of case) it will be removed from the
returned string.

__EXAMPLE__

     mHelp = Informer.gethelptext( 'Informer.gethelptext' ) 

% diplays:
%
% "GETHELPTEXT Return help-text of script|function|method|property as string-vector"
%  ""
% "mHelp = GETHELPTEXT( name )"

     Informer.extracthelpheader( mHelp, "GETHELPTEXT" )

% diplays:
%
% "Return help-text of script|function|method|property as string-vector"

__ETC__

See also
Informer.gethelptext


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : mHelpHeader
- DefiningClass : Informer

---


### extracthelpbody

**Synopsis**: _Return the body of help-text, trimmed of the leading line_ 

      
      mHelpBody = EXTRACTHELPBODY( mHelp ) ;

EXTRACTHELPBODY returns a body of help-text (i.e. from Informer.gethelptext)
and trims it of its leading line of text. (If the text consists solely of the
leading line, then this line is returned.)

__EXAMPLE__
```
% To display the current section of text, without the title line:
mHelpBody = Informer.extracthelpbody( Informer.gethelptext( 'Informer.extracthelpbody' ) )
```

__ETC__
See also

-Informer.gethelptext
-Informer.extracthelpheader


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : mHelpBody
- DefiningClass : Informer

---


### metainfo

**Synopsis**: _EXTRACTMETAINFO Return info (struct) derived from a meta.object_ 

     
     Att = metainfo( MetaObj )

Useful properties of meta[.class/.property/.method] object 'MetaObj' are
copied to the fields of struct Att

TODO : restructuring of meta.Validation.Size


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>true</td></tr>
</table>

- Access : public
- InputNames : Mc
- OutputNames : Att
- DefiningClass : Informer

---


### empty

**Synopsis**: _Returns an empty object array of the given size_ 


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>true</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : E
- DefiningClass : Informer
