# Documentor

**Filetype:** _MATLAB&reg; classdef_

**Synopsis:** _Custom Matlab documentation into markup/down text files_

The DOCUMENTOR class serves to print custom Matlab documentation to file
(e.g. as [Markdown](https://daringfireball.net/projects/markdown/))
while avoiding external dependencies (e.g. [sphinx](https://github.com/sphinx-contrib/matlabdomain))
and tagging syntaxes at odds with Matlab's own style of
[markup](https://www.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html)

Viz., given a list of
[properly](https://www.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html)
commented .m files, a DOCUMENTOR instance outputs simple, readable text files
which are readily hosted online.

### [Example](https://github.com/neuropoly/realtime_shimming/blob/helpDocMd/helpDocMd/doc/Documentor.md)

### Usage

Create a Documentor instance with the list of .m file paths of interest, then
call `printdoc` to write to file:
     
     Dr = Documentor( mFiles ) ;  

     Dr.printdoc ; 

### References

To test how a markdown sample will display when reformatted to HTML:
- [https://daringfireball.net/projects/markdown/dingus](https://daringfireball.net/projects/markdown/dingus)

Re: hosting documentation online:
- [MkDocs](https://www.mkdocs.org/),
- [Github](https://pages.github.com/),
- [ReadTheDocs](https://docs.readthedocs.io/en/stable/)

    Documentation for Documentor
       doc Documentor


### Class Attributes


[table](table)
<table border=1><tr><th>Hidden</th><th>Sealed</th><th>Abstract</th><th>Enumeration</th><th>ConstructOnLoad</th><th>HandleCompatible</th><th>RestrictsSubclassing</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td>[/tr](/tr)
[/table](/table)

- InferiorClasses : [N/A] 
- ContainingPackage : [N/A] 
- EventList : [N/A] 
- EnumerationMemberList : [N/A] 
- Superclasses: handle

- - -
## Properties


### Ext

**Synopsis:** _Default file extensions_

  Default file extensions

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>false</td><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : none
- GetMethod : 
- SetMethod : 
- DefaultValue : [N/A] 
- Validation : [N/A] 
- DefiningClass : Documentor

### Info

**Synopsis:** _Informer instance: Provides the info-content to document a given .m file_

  Informer instance: Provides the info-content to document a given .m file

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : 
- DefaultValue : 
- Validation: 
Class: Informer
Validator functions: 
- DefiningClass : Documentor

### mFiles

**Synopsis:** _List of .m files to document (string scalar or vector of full file paths)_

  List of .m files to document (string scalar or vector of full file paths)

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : Documentor.set.mFiles
- DefaultValue : /Users/ryan/Projects/General/scripts/shim/helpDocMd/src/@Documentor/Documentor.m
- Validation: 
Validator functions: mustBeFile
- DefiningClass : Documentor

### iM

**Synopsis:** _Index of next .m file in mFiles list to document_

  Index of next .m file in mFiles list to document

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : Documentor.set.iM
- DefaultValue : [N/A] 
- Validation: 
Class: uint64
Validator functions: mustBePositive,mustBeInteger
- DefiningClass : Documentor

### isOverwriting

**Synopsis:** _Toggle whether to overwrite existing documentation files_

  Toggle whether to overwrite existing documentation files

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th><th>DefaultValue</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : 
- Validation: 
Validator functions: mustBeNumericOrLogical
- DefiningClass : Documentor

### isSearchRecursive

**Synopsis:** _Toggle whether subdirectories are included in file search (multiple input case only)_

  Toggle whether subdirectories are included in file search (multiple input case only)

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th><th>DefaultValue</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : 
- Validation: 
Validator functions: mustBeNumericOrLogical
- DefiningClass : Documentor

### isSaveRecursive

**Synopsis:** _Recreates original directory tree in dirOut (multiple doc output case only)_

See also HelpDocMd.isSeachRecursive
TODO use mapdirectorytree.m or something to figure out the subdirectory structure to use and change the default to TRUE.
(for now, just dumping all documentation into single folder - dirOutTop

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th><th>DefaultValue</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td><td>false</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : 
- Validation: 
Validator functions: mustBeNumericOrLogical
- DefiningClass : Documentor

### dirOutTop

**Synopsis:** _Output parent directory for the doc files_


See also HelpDocMd.isSaveRecursive

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : Documentor.get.dirOutTop
- SetMethod : Documentor.set.dirOutTop
- DefaultValue : 
- Validation: 
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

### mdDoc

**Synopsis:** _Reformated documentation_

  Reformated documentation

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : Documentor.get.mdDoc
- SetMethod : 
- DefaultValue : 
- Validation: 
Class: string
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

### syntax

**Synopsis:** _String specifier for output syntax: "mkd" (for Mkdocs markdown), "mat" (for MATLAB markup)_

The sole difference between "mkd" and "mat" (for now) is that "mkd" will
reformat the style of any embedded links in the comments.

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : 
- DefaultValue : mkd
- Validation: 
Class: string
Validator functions: @(syntax)mustBeMember(syntax,["mat","mkd"])
- DefiningClass : Documentor

### isDetailed

**Synopsis:** _Toggles between basic/user (=false) and detailed/developer documentation (=true) [default: true]_

When false, classes and class members with private, protected, or hidden
attributes are excluded from the output documentation. [default = true]

TODO: implementation!

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th><th>DefaultValue</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : 
- Validation: 
Validator functions: mustBeBoolean
- DefiningClass : Documentor

### extOut

**Synopsis:** _Output file extension (default = ".md")_

  Output file extension (default = ".md")

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : public
- SetAccess : public
- GetMethod : 
- SetMethod : 
- DefaultValue : .md
- Validation: 
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

### mDir

**Synopsis:** _parent folder of mFiles(iM)_

  parent folder of mFiles(iM)  

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : private
- SetAccess : private
- GetMethod : 
- SetMethod : 
- DefaultValue : /Users/ryan/Projects/General/scripts/shim/helpDocMd/src/@Documentor
- Validation: 
Validator functions: mustBeFolder
- DefiningClass : Documentor

### dirInTop

**Synopsis:** _top directory of src mFiles_

  top directory of src mFiles

[table](table)
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th>[/tr](/tr)
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>true</td>[/tr](/tr)
[/table](/table)

- GetAccess : private
- SetAccess : private
- GetMethod : Documentor.get.dirInTop
- SetMethod : 
- DefaultValue : 
- Validation: 
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

---
## Methods


---


### Documentor

**Synopsis**: _Custom Matlab documentation into markup/down text files_ 

The DOCUMENTOR class serves to print custom Matlab documentation to file
(e.g. as [Markdown](https://daringfireball.net/projects/markdown/))
while avoiding external dependencies (e.g. [sphinx](https://github.com/sphinx-contrib/matlabdomain))
and tagging syntaxes at odds with Matlab's own style of
[markup](https://www.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html)

Viz., given a list of
[properly](https://www.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html)
commented .m files, a DOCUMENTOR instance outputs simple, readable text files
which are readily hosted online.

### [Example](https://github.com/neuropoly/realtime_shimming/blob/helpDocMd/helpDocMd/doc/Documentor.md)

### Usage

Create a Documentor instance with the list of .m file paths of interest, then
call `printdoc` to write to file:
     
     Dr = Documentor( mFiles ) ;  

     Dr.printdoc ; 

### References

To test how a markdown sample will display when reformatted to HTML:
- [https://daringfireball.net/projects/markdown/dingus](https://daringfireball.net/projects/markdown/dingus)

Re: hosting documentation online:
- [MkDocs](https://www.mkdocs.org/),
- [Github](https://pages.github.com/),
- [ReadTheDocs](https://docs.readthedocs.io/en/stable/)


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : pathIn, Params
- OutputNames : Dr
- DefiningClass : Documentor

---


### printdoc

**Synopsis**: _Write documentation to file_ 

### Syntax
     
     [filename] = PRINTDOC( Self ) 

Writes the contents of `Self.mdDoc` to `filename`.
To overwrite an existing file, first set `Self.isOverwriting = true`


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : Dr
- OutputNames : pathOut
- DefiningClass : Documentor

---


### findfilestodocument

**Synopsis**: _Return list of .m files to document from directory search_ 

FINDFILESTODOCUMENT searches a directory for .m files and then removes any
class method files from the list (methods are included as part of the
overall class documentation).

### Syntax
     
     mFiles = Dr.FINDFILESTODOCUMENT( pathIn )

### Implementation details

FINDFILESTODOCUMENT wraps to findfiles with the function call:
      
     [~,mFiles] = findfiles( folder, "*.m", Dr.isSearchRecursive ) ;

.m file types are then determined using Informer.mfiletype( mFiles ) and, if
present in the list, methods .m files are removed.


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : Dr, pathIn
- OutputNames : mFiles
- DefiningClass : Documentor

---


### tableattributes

**Synopsis**: _Return html-table of class/classmember attributes_ 


tableStr = GETHELPTEXT( Attributes )


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr, Attributes
- OutputNames : tableStr
- DefiningClass : Documentor

---


### documentclassproperties

**Synopsis**: _Return string vector of class property documentation_ 

 DOCUMENTCLASSPROPERTIES Return string vector of class property documentation


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

---


### documentclassmethods

**Synopsis**: _Return string vector of class method documentation_ 

 DOCUMENTCLASSMETHODS Return string vector of class method documentation


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

---


### documentclassdef

**Synopsis**: _Return string vector of class documentation_ 

DOCUMENTCLASSDEF documents basic class attributes followed by class member
documentation (courtesy of calls to Documentor.documentclassproperties and
Documentator.documentclassmethods)


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

---


### documentbasic

**Synopsis**: _Return string vector of rudimentary documentation_ 

### Syntax
     
     [docStr] = DOCUMENTBASIC( Self )
     [docStr] = DOCUMENTBASIC( Self, Att )
     [docStr] = DOCUMENTBASIC( Self, Att, headingLevel )

When called without a second argument, DOCUMENTBASIC derives the following details
from `Self.Info.Attributes` to return the documentation string vector `docStr`:
- Name : of the script, function, or class
- mType: script, function, or class file type
- Description: header line from the help/documentation
- DetailedDescription: body of the help/documentation

When called with a second argument, DOCUMENTBASIC works similarly to the
above case, however, details are instead derived from attributes-struct `Att`
(mType may be omitted in this case, but the other field names must be
present.)

Optional 3rd argument is a scalar integer (= 0,1,2,3,4,5, or 6) [default=1]
indicating the number of '#' signs to precede the 'Name' value (docStr's
first element). for markdown syntax, this controls the heading level.


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : varargin
- OutputNames : docStr
- DefiningClass : Documentor

---


### markuptodown

**Synopsis**: _Replace MATLAB Markup elements with corresponding Mkdocs-style Markdown_ 

Reformat links and link-text to Markdown syntax:
i.e. MATLAB markup uses: [texttodisplay](https://www.thisSite.com)
whereas Markdown uses: [text to display](https://www.thisSite.com)

### Syntax ###

[ mdDocStr ] = MARKUPTODOWN( muDocStr )

### NOTE/TODO ###

The implementation is simplistic: basically works for weblinks, but needs
to be elaborated for local links to custom functions & classes: either tags
or relative paths could be used for Mkdocs build.

Moreover, it doesn't distinguish between links and embedded HTML, so it
messes up the latter (see: Documentor.tableattributes)


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : muDocStr
- OutputNames : mdDocStr
- DefiningClass : Documentor

---


### documentfunction

**Synopsis**: _adds function-specific info to documentation_ 

(NOTE: for now, this is just nArgin/nArgout but this should be elaborated
in Informer.m -- e.g. by parsing the function arguments block when it exists)


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : varargin
- OutputNames : varargout
- DefiningClass : Documentor

---


### empty

**Synopsis**: _Returns an empty object array of the given size_ 


#### Attributes:

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>true</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : varargin
- OutputNames : E
- DefiningClass : Documentor

---

### eq
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### ne
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### lt
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### gt
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### le
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### ge
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### delete
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### isvalid
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### findprop
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### notify
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### notify
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### addlistener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### listener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### addlistener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### listener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### addlistener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### listener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### addlistener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### listener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### addlistener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### listener
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]

---

### findobj
[ _built-in method derived from class_ **handle** ]
For more info, see MATLAB documentation]
