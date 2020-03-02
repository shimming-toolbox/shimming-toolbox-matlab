# Documentor # (a MATLAB(R) classdef)

_Custom MATLAB documentation into markup/down text files_

### Description ###



Writes *thorough* Matlab documentation as simple, readable
<https://daringfireball.net/projects/markdown/ Markdown> text which is

1. Readily hosted online (e.g. <https://www.mkdocs.org/ MkDocs>,
<https://pages.github.com/>, <https://docs.readthedocs.io/en/stable/>)

2. Does not require additional dependencies or different syntax/tagging
from Matlab's own markup style (e.g. sphinx)

### Basic Usage ###

1. User creates a Documentor instance with the list of .m file paths to
document:

     Dr = Documentor( mFiles ) ;

2. To create the .md documentation, the user calls:

     Dr.write ; 

### Example Output ###

(TODO: add output to github page or readthedocs)

To see the final output online, see <https://ADD_URL.com this> %

### References ###

To test how the .md output will appear once reformatted to HTML:
<https://daringfireball.net/projects/markdown/dingus>

    Documentation for Documentor
       doc Documentor




### Class Attributes ###


<table>
<table border=1><tr><th>Hidden</th><th>Sealed</th><th>Abstract</th><th>Enumeration</th><th>ConstructOnLoad</th><th>HandleCompatible</th><th>RestrictsSubclassing</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>false</td></tr>
</table>

- InferiorClasses : [N/A] 
- ContainingPackage : [N/A] 
- EventList : [N/A] 
- EnumerationMemberList : [N/A] 
- Superclasses: handle
*dirInTop*

_top directory of src mFiles_



<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>true</td></tr>
</table>

- GetAccess : private
- SetAccess : private
- GetMethod : Documentor.get.dirInTop
- SetMethod : 
- DefaultValue : 
- Validation: 
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor
---
### Methods ###


#### Documentor #### : _Custom MATLAB documentation into markup/down text files_
Description: 

Writes *thorough* Matlab documentation as simple, readable
<https://daringfireball.net/projects/markdown/ Markdown> text which is

1. Readily hosted online (e.g. <https://www.mkdocs.org/ MkDocs>,
<https://pages.github.com/>, <https://docs.readthedocs.io/en/stable/>)

2. Does not require additional dependencies or different syntax/tagging
from Matlab's own markup style (e.g. sphinx)

### Basic Usage ###

1. User creates a Documentor instance with the list of .m file paths to
document:

     Dr = Documentor( mFiles ) ;

2. To create the .md documentation, the user calls:

     Dr.write ; 

### Example Output ###

(TODO: add output to github page or readthedocs)

To see the final output online, see <https://ADD_URL.com this> %

### References ###

To test how the .md output will appear once reformatted to HTML:
<https://daringfireball.net/projects/markdown/dingus>

- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : pathIn, Params
- OutputNames : Dr
- DefiningClass : Documentor

#### write #### : _Write documentation to file_
Description: 

WRITE creates (or, optionally, overwrites) a .md file and writes to it the
contents of Self.mdDoc.

### Syntax ###

[pathOut] = write( Self )

To enable overwriting of existing files, set Self.isOverwriting = true

- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : Dr
- OutputNames : pathOut
- DefiningClass : Documentor

#### findfilestodocument #### : _Return list of .m files to document from directory search_
Description: 

FINDFILESTODOCUMENT searches a directory for .m files and then removes any
class method files from the list (methods are documented as part of the
overall class documentation).

### Syntax ###

mFiles = Dr.FINDFILESTODOCUMENT( pathIn )

### Implementation Details ###

FINDFILESTODOCUMENT wraps to findfiles with the function call:
    
    [~,mFiles] = findfiles( folder, "*.m", Dr.isSearchRecursive ) ;

.m file types are then determined using Informer.mfiletype( mFiles ) and, if
present in the list, methods .m files are removed.
     
TODO give user option of documenting retaining these methods files if for
some reason there is a need to document them independently of the class?

- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : Dr, pathIn
- OutputNames : mFiles
- DefiningClass : Documentor

#### documentclassproperties #### : _Return string vector of class property documentation_
Description: 

- Access : private
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

#### documentclassmethods #### : _Return string vector of class method documentation_
Description: 

- Access : private
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

#### documentclassdef #### : _Return string vector of class documentation_
Description: 

DOCUMENTCLASSDEF documents basic class attributes followed by class member
documentation (courtesy of calls to Documentor.documentclassproperties and
Documentator.documementclassmethods)

- Access : private
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

#### documentbasic #### : _Return string vector of rudimentary documentation for all .m file types_
Description: 

DOCUMENTBASIC Documents the following .m file details:
- Name : of the script, function, or class
- Type: script, function, or class
- Description: header line from the help/documentation
- DetailedDescription: body of the help/documentation

- Access : private
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : Dr
- OutputNames : docStr
- DefiningClass : Documentor

#### markuptodown #### : _Replace MATLAB Markup elements with corresponding Mkdocs-style Markdown_
Description: 

Reformat links and link-text to Markdown syntax:
i.e. MATLAB markup uses: <https://www.thisSite.com text to display>
whereas Markdown uses: [text to display](https://www.thisSite.com)

### Syntax ###

[ mdDocStr ] = MARKUPTODOWN( muDocStr )

### NOTE/TODO ###

The implementation works for weblinks, but will need to be elaborated
for local links to custom functions & classes: either tags or relative
paths could be used for Mkdocs build.

- Access : public
- Static : true
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : muDocStr
- OutputNames : mdDocStr
- DefiningClass : Documentor

#### documentfunction #### : _adds function-specific info to documentation_
Description: 
(NOTE: for now, this is just nArgin/nArgout but this should be elaborated
in Informer.m -- e.g. by parsing the function arguments block when it exists)

- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : varargin
- OutputNames : varargout
- DefiningClass : Documentor

#### empty #### : _Returns an empty object array of the given size_
- Access : public
- Static : true
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : true
- InputNames : varargin
- OutputNames : E
- DefiningClass : Documentor

#### eq ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### ne ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### lt ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### gt ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### le ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### ge ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### delete ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### isvalid ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### findprop ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### notify ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### notify ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### addlistener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### listener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### addlistener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### listener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### addlistener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### listener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### addlistener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### listener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### addlistener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### listener ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]

#### findobj ####[ _built-in method derived from *handle* class_ ]
For more info, see MATLAB documentation]
