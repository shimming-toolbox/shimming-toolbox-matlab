# Documentor

**Filetype:** _MATLAB&reg; classdef_

**Synopsis:** _Custom MATLAB documentation into markup/down text files_

The DOCUMENTOR class serves to write custom Matlab documentation as simple,
readable text files (i.e. [Markdown](https://daringfireball.net/projects/markdown/)).
Specifically, it serves to produce documentation that is readily hosted
online (e.g. <https://www.mkdocs.org/ MkDocs>, [https://pages.github.com/](https://pages.github.com/),
[https://docs.readthedocs.io/en/stable/](https://docs.readthedocs.io/en/stable/)) without the need for external
dependencies or tagging syntaxes that differ from Matlab's own markup style
(e.g. [sphinx](https://github.com/sphinx-contrib/matlabdomain)).

### Usage

Create a Documentor instance with the list of .m file paths of interest, then
call `write` to write to file:

`
Dr = Documentor( mFiles ) ;
Dr.write ;
`

### Example

[helpDocMd](https://github.com/neuropoly/realtime_shimming/blob/helpDocMd/helpDocMd/doc/Documentor.md)

(TODO: replace with link to github page or readthedocs page with the complete
documentation for the entire helpDocMd code library)

### References

To test how the .md output will appear once reformatted to HTML:
[https://daringfireball.net/projects/markdown/dingus](https://daringfireball.net/projects/markdown/dingus)

    Documentation for Documentor
       doc Documentor


### Class Attributes ###


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
## Properties ##


#### Ext ####

_Default file extensions_

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

#### Info ####

_Informer instance: Provides the info-content to document a given .m file_

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

#### mFiles ####

_List of .m files to document (string scalar or vector of full file paths)_

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

#### iM ####

_Index of next .m file in mFiles list to document_

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

#### isOverwriting ####

_Toggle whether to overwrite existing documentation files_

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

#### isSearchRecursive ####

_Toggle whether subdirectories are included in file search (multiple input case only)_

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

#### isSaveRecursive ####

_Recreates original directory tree in dirOut (multiple doc output case only)_

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

#### dirOutTop ####

_Output parent directory for the doc files_


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

#### mdDoc ####

_Reformated documentation_

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

#### syntax ####

_String specifier for output syntax: "mkd" (for Mkdocs markdown), "mat" (for MATLAB markup)_

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

#### isDetailed ####

_Toggles between basic/user (=false) and detailed/developer documentation (=true) [default: true]_

When false, classes and class members with private, protected, or hidden
attributes are excluded from the output documentation. [default = true]

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

#### extOut ####

_Output file extension (default = ".md")_

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

#### mDir ####

_parent folder of mFiles(iM)_

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

#### dirInTop ####

_top directory of src mFiles_

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
## Methods ##


---


### Documentor ###


 _Custom MATLAB documentation into markup/down text files_ 

Description: 
The DOCUMENTOR class serves to write custom Matlab documentation as simple,
readable text files (i.e. [Markdown](https://daringfireball.net/projects/markdown/)).
Specifically, it serves to produce documentation that is readily hosted
online (e.g. <https://www.mkdocs.org/ MkDocs>, [https://pages.github.com/](https://pages.github.com/),
[https://docs.readthedocs.io/en/stable/](https://docs.readthedocs.io/en/stable/)) without the need for external
dependencies or tagging syntaxes that differ from Matlab's own markup style
(e.g. [sphinx](https://github.com/sphinx-contrib/matlabdomain)).

### Usage

Create a Documentor instance with the list of .m file paths of interest, then
call `write` to write to file:

`
Dr = Documentor( mFiles ) ;
Dr.write ;
`

### Example

[helpDocMd](https://github.com/neuropoly/realtime_shimming/blob/helpDocMd/helpDocMd/doc/Documentor.md)

(TODO: replace with link to github page or readthedocs page with the complete
documentation for the entire helpDocMd code library)

### References

To test how the .md output will appear once reformatted to HTML:
[https://daringfireball.net/projects/markdown/dingus](https://daringfireball.net/projects/markdown/dingus)


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : pathIn, Params
- OutputNames : Dr
- DefiningClass : Documentor

---


### write ###


 _Write documentation to file_ 

Description: 
### Syntax
     
     [pathOut] = write( Self ) 

writes the contents of Self.mdDoc to file. To _overwrite_ an existing file,
before calling WRITE, set
      Self.isOverwriting = true


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : Dr
- OutputNames : pathOut
- DefiningClass : Documentor

---


### findfilestodocument ###


 _Return list of .m files to document from directory search_ 

Description: 
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


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : Dr, pathIn
- OutputNames : mFiles
- DefiningClass : Documentor

---


### tableattributes ###


 _Return html-table of class/classmember attributes_ 

Description: 

tableStr = GETHELPTEXT( Attributes )


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr, Attributes
- OutputNames : tableStr
- DefiningClass : Documentor

---


### documentclassproperties ###


 _Return string vector of class property documentation_ 

Description: 
 DOCUMENTCLASSPROPERTIES Return string vector of class property documentation


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

---


### documentclassmethods ###


 _Return string vector of class method documentation_ 

Description: 
 DOCUMENTCLASSMETHODS Return string vector of class method documentation


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

---


### documentclassdef ###


 _Return string vector of class documentation_ 

Description: 
DOCUMENTCLASSDEF documents basic class attributes followed by class member
documentation (courtesy of calls to Documentor.documentclassproperties and
Documentator.documementclassmethods)


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

---


### documentbasic ###


 _Return string vector of rudimentary documentation for all .m file types_ 

Description: 
DOCUMENTBASIC Documents the following .m file details:
- Name : of the script, function, or class
- Type: script, function, or class
- Description: header line from the help/documentation
- DetailedDescription: body of the help/documentation


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : private
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

---


### markuptodown ###


 _Replace MATLAB Markup elements with corresponding Mkdocs-style Markdown_ 

Description: 
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


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : muDocStr
- OutputNames : mdDocStr
- DefiningClass : Documentor

---


### documentfunction ###


 _adds function-specific info to documentation_ 

Description: 
(NOTE: for now, this is just nArgin/nArgout but this should be elaborated
in Informer.m -- e.g. by parsing the function arguments block when it exists)


#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : varargin
- OutputNames : varargout
- DefiningClass : Documentor

---


### empty ###


 _Returns an empty object array of the given size_ 



#### Attributes: ####

[table](table)
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th>[/tr](/tr)
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>true</td>[/tr](/tr)
[/table](/table)

- Access : public
- InputNames : varargin
- OutputNames : E
- DefiningClass : Documentor

---

### eq ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### ne ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### lt ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### gt ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### le ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### ge ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### delete ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### isvalid ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### findprop ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### notify ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### notify ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### addlistener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### listener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### addlistener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### listener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### addlistener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### listener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### addlistener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### listener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### addlistener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### listener ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

---

### findobj ###
[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]
