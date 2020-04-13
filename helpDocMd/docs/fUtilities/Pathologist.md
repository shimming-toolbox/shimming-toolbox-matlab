# Pathologist

**Filetype:** _MATLAB&reg; classdef_

**Synopsis:** _File path utility class_

      
      PathO = Pathologist( pathIn ) ;

Recasts the array of file paths `pathIn` as a `Pathologist` object,
containing simple+convenient methods for handling paths to files and folders.

Since a path designation in MATLAB can take one of three forms (char-,
string-, or cellstr-array), a primary purpose the class is merely to provide
a convenient means of handling the different implementations. Namely, a set
of input paths will always, first, be typecast as a string and stored as the
property `data`. After processing via any `Pathologist` methods, paths are
returned in the form of their input (see Pathologist.returnasinput).

**NOTE** The name `pathtool` already being taken by a builtin MATLAB
function, in keeping with the `Doc. Md.` theme and the unfortunately current
zeitgeist due to COVID-19, "Pathologist" seems an appropriate moniker.

    Documentation for Pathologist
       doc Pathologist


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


### typeIn

**Synopsis:** _Input data-type (aka 'class')_

  Input data-type (aka 'class') 

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td></tr>
</table>

- GetAccess : private
- SetAccess : private
- PartialMatchPriority : [N/A] 
- GetMethod : 
- SetMethod : 
- DefaultValue : char
- Validation: 
Class: string
Validator functions: @(typeIn)mustBeMember(typeIn,["string","char","cell"])
- DefiningClass : Pathologist

### sizeIn

**Synopsis:** _Pathologist/sizeIn is a property._

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td></tr>
</table>

- GetAccess : private
- SetAccess : private
- PartialMatchPriority : [N/A] 
- GetMethod : 
- SetMethod : 
- DefaultValue : [N/A] 
- Validation : [N/A] 
- DefiningClass : Pathologist

### data

**Synopsis:** _Path data: real or fictitious files or folders_

  Path data: real or fictitious files or folders

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : 
- SetMethod : 
- DefaultValue : ./
- Validation: 
Validator functions: mustBeString
- DefiningClass : Pathologist

### typeOut

**Synopsis:** _Sets the return datatype (class) of path info as "string" "char" or "cell"_


(Default is whatever data type was used to initialize the object.)

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : 
- SetMethod : 
- DefaultValue : char
- Validation: 
Class: string
Validator functions: @(typeOut)mustBeMember(typeOut,["string","char","cell"])
- DefiningClass : Pathologist

### files

**Synopsis:** _list of files among Path.data_

  list of files among Path.data

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : Pathologist.get.files
- SetMethod : 
- DefaultValue : 
- Validation : [N/A] 
- DefiningClass : Pathologist

### folders

**Synopsis:** _list of folders among Path.data_

  list of folders among Path.data

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : Pathologist.get.folders
- SetMethod : 
- DefaultValue : 
- Validation : [N/A] 
- DefiningClass : Pathologist

### subfolders

**Synopsis:** _of each folder element of Path.data_

  subfolders of each folder element of Path.data

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : Pathologist.get.subfolders
- SetMethod : 
- DefaultValue : 
- Validation : [N/A] 
- DefiningClass : Pathologist

### nPaths

**Synopsis:** _Number of input paths_

  Number of input paths

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : Pathologist.get.nPaths
- SetMethod : 
- DefaultValue : 
- Validation : [N/A] 
- DefiningClass : Pathologist

### parentDir

**Synopsis:** _Path(s) to immediate parent folder(s) (single level)_

Main doc. entry (for now?): parent.m

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : Pathologist.get.parentDir
- SetMethod : 
- DefaultValue : 
- Validation : [N/A] 
- DefiningClass : Pathologist

### baseDir

**Synopsis:** _Deepest common parent folder (the "fork")_

`if Path.nPaths == 1, Path.baseDir = Path.parentDir`
Main doc. entry (for now?): parent.m

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td></tr>
</table>

- GetAccess : public
- SetAccess : public
- PartialMatchPriority : [N/A] 
- GetMethod : Pathologist.get.baseDir
- SetMethod : 
- DefaultValue : 
- Validation: 
Validator functions: 
- DefiningClass : Pathologist

---
## Methods


---


### isvalid

**Synopsis**: _Returns true for elements of Path.data that refer to files or folders_ 


`Path.isValid = [ Path.isfile | Path.isfolder ]`


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : Path
- OutputNames : isValid
- DefiningClass : Pathologist

---


### Pathologist

**Synopsis**: _File path utility class_ 

      
      PathO = Pathologist( pathIn ) ;

Recasts the array of file paths `pathIn` as a `Pathologist` object,
containing simple+convenient methods for handling paths to files and folders.

Since a path designation in MATLAB can take one of three forms (char-,
string-, or cellstr-array), a primary purpose the class is merely to provide
a convenient means of handling the different implementations. Namely, a set
of input paths will always, first, be typecast as a string and stored as the
property `data`. After processing via any `Pathologist` methods, paths are
returned in the form of their input (see Pathologist.returnasinput).

**NOTE** The name `pathtool` already being taken by a builtin MATLAB
function, in keeping with the `Doc. Md.` theme and the unfortunately current
zeitgeist due to COVID-19, "Pathologist" seems an appropriate moniker.


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : pathIn
- OutputNames : Path
- DefiningClass : Pathologist

---


### fileparts

**Synopsis**: _Overloaded MATLAB function: [filepath,name,ext] = arrayfun( @fileparts, Path.data ) ;_ 

 FILEPARTS Overloaded MATLAB function: [filepath,name,ext] = arrayfun( @fileparts, Path.data ) ;


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : Path
- OutputNames : folder, name, ext
- DefiningClass : Pathologist

---


### exist

**Synopsis**: _Overloaded MATLAB function, returnArg = arrayfun( @exist, Path.data ) ;_ 

 EXIST Overloaded MATLAB function, returnArg = arrayfun( @exist, Path.data ) ;


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : Path
- OutputNames : mExist
- DefiningClass : Pathologist

---


### isfolder

**Synopsis**: _Overloaded MATLAB function, returnArg = isfolder( Path.data ) ;_ 

 ISFOLDER Overloaded MATLAB function, returnArg = isfolder( Path.data ) ;


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : Path
- OutputNames : isDir
- DefiningClass : Pathologist

---


### isfile

**Synopsis**: _ISFOLDER Overloaded MATLAB function, returnArg = isfile( Path.data ) ;_ 

 ISFOLDER Overloaded MATLAB function, returnArg = isfile( Path.data ) ;


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : Path
- OutputNames : isFilepath
- DefiningClass : Pathologist

---


### subsref

**Synopsis**: _Overloaded MATLAB function_ 

 SUBSREF Overloaded MATLAB function


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td></tr>
</table>

- Access : public
- InputNames : Path, S
- OutputNames : pathOut
- DefiningClass : Pathologist

---


### returnasinput

**Synopsis**: _Resize and recast a path value to correspond with user input_ 

     
      [pathOut] = returnasinput( Path )
      [pathOut] = returnasinput( Path, pathIn )
      [pathOut1, pathOut2, ...] = returnasinput( Path, pathIn1, pathIn2,... )

Reshapes `pathIn` according to `Path.sizeIn` when the 2 have the same number
of elements and recasts it as `Path.typeIn` to return `pathOut`.

When called with the single argument, `pathIn` is assigned the value of
`Path.data`.

When called with > 2 inputs, resizing + recasting is applied successively to
define the respective returns.


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : private
- InputNames : Path, varargin
- OutputNames : varargout
- DefiningClass : Pathologist

---


### whatr

**Synopsis**: _Recursive` what()`_ 

      
      Info = whatr( baseDir, isExcludingHidden )

Calls MATLAB function [what] for `baseDir` and each of its subdirectories,
and arrays the results to return `Info`.

__INPUTS__
      
    `baseDir=["."]`
      The parent/top directory path as a string scalar or character vector.

    `isExcludingHidden=[true|1]`
      A Boolean toggle: when set to false, hidden folders are included in the returned list.

__ETC__

To help remember the function name, think "what -r" on the commandline, and
"what files are there?" in English.

- [what](https://www.mathworks.com/help/matlab/ref/what.html)

See also
WHAT


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : Info
- DefiningClass : Pathologist

---


### mapdirectorytree

**Synopsis**: _Returns list of subdirectories_ 

     
      dirs = MAPDIRECTORYTREE( baseDir, isReturnRelative, isExcludingHidden )

Returns paths to the subdirectories of `baseDir` as elements of the string
column vector `dirs`.

__INPUTS__

    baseDir=["."]
      Parent/top directory as string scalar

    isReturnRelative=[true|1]
      Boolean toggle: Returns the subdirectory paths as relative to `baseDir` 
      when true, or in absolute terms when false.

    isExcludingHidden=[true|1]
      Boolean toggle includes hidden folders in the returned list when set to 0
      (false).

__ETC__

See also
FILEATTRIB, DIR


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : dirs
- DefiningClass : Pathologist

---


### typecast

**Synopsis**: _Convert between string, char, cellstr_ 

      [txtOut] = typecast( txtIn ) 
      [txtOut] = typecast( txtIn, castAs ) 

Converts an array of text `txtIn` to the class specified in `castAs`, which
may be `"string"`, `"char"` or `"cell"` (or, equivalently, `"cellstr"`).
When called with a single input, the function returns the text as a string.
As with [cellstr], trailing white space is deleted when converting to
string from char.

__ETC__

Matlab has three options for text-containers:

1. `char`, which is the most primitive and least flexible.

2. `cell`, which is a general container for any data type.

3. `string` which has the virtue of being much more flexible than `char`
while still being specialized for text, unlike `cell`. However, as the most
recently introduced (2016b) of the three, it is often necessary to convert
to one of the other two.

[string](https://www.mathworks.com/help/matlab/ref/string.html)
[cellstr](https://www.mathworks.com/help/matlab/ref/cellstr.html)

See also
CHAR,CELL,CELLSTR,STRING


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : txtOut
- DefiningClass : Pathologist

---


### abs

**Synopsis**: _Validate and convert file system paths_ 

     
     [pathOut, pathType] = abs( pathIn )  

Wraps to the Matlab function `fileattrib` to check each element of the input
array `pathIn` for valid file system paths (which can be relative or
abs/absolute). It returns two arrays of the same type (string, char, or
cellstr) and size (arbitrary) as the input.

The possibilities for a given return element `(i)` are:

`pathIn(i)`  | `pathOut(i)`| `pathType(i)`
-------------|-------------|-------------------------------------
is a file   | abs path   | "file"
is a folder | abs path   | "directory"
else        |   ""        | the error message from [fileattrib]

**Note** for invalid input paths, the error message returned from
`fileattrib( pathIn(i) )` will presumably be 'No such file or directory.'
However the MATLAB documentation is unclear whether other possibilities exist
for the case where the function is called with a single argument, as is done
here. For more info, refer to the MATLAB documentation:
[fleattrib](https://www.mathworks.com/help/matlab/ref/fileattrib.html/)

See also
FILEATTRIB


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : varargin
- OutputNames : pathOut, pathType
- DefiningClass : Pathologist

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
- DefiningClass : Pathologist
