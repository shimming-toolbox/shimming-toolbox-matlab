# Info.Name # (a Matlab classdef)

_Custom MATLAB documentation into markup/down text files_

### Description ###


Writes *thorough* Matlab documentation as simple, readable Markdown text
<https://daringfireball.net/projects/markdown/> which is

1. readily hosted online (e.g. <https://www.mkdocs.org>,
<https://pages.github.com/>, <https://docs.readthedocs.io/en/stable/>)

2. does not require additional dependencies or different syntax/tagging
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

### Attributes ###
- Hidden : false
- Sealed : false
- Abstract : false
- Enumeration : false
- ConstructOnLoad : false
- HandleCompatible : true
- InferiorClasses : [N/A]
- ContainingPackage : [N/A]
- RestrictsSubclassing : false
- EventList : [N/A]
- EnumerationMemberList : [N/A]
- Superclasses: handle

### Properties ###


*Ext* : _Default file extensions_
- GetAccess : public
- SetAccess : none
- Dependent : false
- Constant : true
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : false
- NonCopyable : true
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : [N/A]
- Validation : [N/A]
- DefiningClass : Documentor

*Info* : _Informer object instance: Provides the information content to document a_
Description:
given .m file
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod :
- HasDefault : false
- DefaultValue :
- Validation:
Class: Informer
Validator functions:
- DefiningClass : Documentor

*mFiles* : _List of .m files to document (string scalar or vector of full file paths)_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod : Documentor.set.mFiles
- HasDefault : true
- DefaultValue : /Users/ryan/Projects/General/scripts/shim/helpDocMd/src/@Documentor/Documentor.m
- Validation:
Validator functions: mustBeFile
- DefiningClass : Documentor

*iM* : _Index of next .m file in mFiles list to document_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod : Documentor.set.iM
- HasDefault : true
- DefaultValue : [N/A]
- Validation:
Class: uint64
Validator functions: mustBePositive,mustBeInteger
- DefiningClass : Documentor

*isOverwriting* : _Toggle whether to overwrite existing documentation files_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : false
- Validation:
Validator functions: mustBeNumericOrLogical
- DefiningClass : Documentor

*isSearchRecursive* : _Toggle whether subdirectories are included in file search (multiple input case only)_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : true
- Validation:
Validator functions: mustBeNumericOrLogical
- DefiningClass : Documentor

*isSaveRecursive* : _Recreates original directory tree in dirOut (multiple doc output case only)_
Description:

See also HelpDocMd.isSeachRecursive
TODO use mapdirectorytree.m or something to figure out the subdirectory structure to use and change the default to TRUE.
(for now, just dumping all documentation into single folder - dirOutTop
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : false
- Validation:
Validator functions: mustBeNumericOrLogical
- DefiningClass : Documentor

*dirOutTop* : _Output parent directory for the doc files_
Description:

See also HelpDocMd.isSaveRecursive
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod : Documentor.get.dirOutTop
- SetMethod : Documentor.set.dirOutTop
- HasDefault : true
- DefaultValue :
- Validation:
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

*mdDoc* : _Reformated documentation_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod : Documentor.get.mdDoc
- SetMethod :
- HasDefault : true
- DefaultValue :
- Validation:
Class: string
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

*syntax* : _String specifier for output syntax: "mkd" (for Mkdocs markdown), "mat" (for MATLAB markup)_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : mkd
- Validation:
Class: string
Validator functions: @(syntax)mustBeMember(syntax,["mat","mkd"])
- DefiningClass : Documentor

*extIn* : _Input/Matlab file extension_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : .m
- Validation:
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

*nameIn* : _Input names (without directory path or file extension)_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : Documentor
- Validation:
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

*extOut* : _Output file extension (default = ".md")_
- GetAccess : public
- SetAccess : public
- Dependent : false
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : true
- NonCopyable : false
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : .md
- Validation:
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

*mDir* : _parent folder of mFiles(iM)_
- GetAccess : private
- SetAccess : private
- Dependent : true
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : false
- NonCopyable : true
- GetMethod :
- SetMethod :
- HasDefault : true
- DefaultValue : /Users/ryan/Projects/General/scripts/shim/helpDocMd/src/@Documentor
- Validation:
Validator functions: mustBeFolder
- DefiningClass : Documentor

*dirInTop* : _top directory of src mFiles_
- GetAccess : private
- SetAccess : private
- Dependent : true
- Constant : false
- Abstract : false
- Transient : false
- Hidden : false
- GetObservable : false
- SetObservable : false
- AbortSet : false
- NonCopyable : true
- GetMethod : Documentor.get.dirInTop
- SetMethod :
- HasDefault : true
- DefaultValue :
- Validation:
Validator functions: mustBeStringOrChar
- DefiningClass : Documentor

### Methods ###


####write#### : _Write doc contents to file_
Description:

[pathOut] = write( Self )
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : Dr
- OutputNames : pathOut
- DefiningClass : Documentor

####docclassmethods#### : _DOCCLASSMETHODS_
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : Dr
- OutputNames : membersDoc
- DefiningClass : Documentor

####Documentor#### : _Custom MATLAB documentation into markup/down text files_
Description:

Writes *thorough* Matlab documentation as simple, readable Markdown text
<https://daringfireball.net/projects/markdown/> which is

1. readily hosted online (e.g. <https://www.mkdocs.org>,
<https://pages.github.com/>, <https://docs.readthedocs.io/en/stable/>)

2. does not require additional dependencies or different syntax/tagging
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

####documentclassattributes#### : _DOCUMENTCLASSATTRIBUTES_
- Access : public
- Static : true
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : varargin
- OutputNames : info
- DefiningClass : Documentor

####findfilestodocument#### : _Return list of .m files to document from directory search_
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

####empty#### : _Returns an empty object array of the given size_
- Access : public
- Static : true
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : true
- InputNames : varargin
- OutputNames : E
- DefiningClass : Documentor

####eq#### : _== (EQ)   Test handle equality._
Description:
Handles are equal if they are handles for the same object.

H1 == H2 performs element-wise comparisons between handle arrays H1 and
H2.  H1 and H2 must be of the same dimensions unless one is a scalar.
The result is a logical array of the same dimensions, where each
element is an element-wise equality result.

If one of H1 or H2 is scalar, scalar expansion is performed and the
result will match the dimensions of the array that is not scalar.

TF = EQ(H1, H2) stores the result in a logical array of the same
dimensions.

See also HANDLE, HANDLE/GE, HANDLE/GT, HANDLE/LE, HANDLE/LT, HANDLE/NE

Documentation for handle/eq
doc handle.eq
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : A, B
- OutputNames : TF
- DefiningClass : handle

####ne#### : _~= (NE)   Not equal relation for handles._
Description:
Handles are equal if they are handles for the same object and are
unequal otherwise.

H1 ~= H2 performs element-wise comparisons between handle arrays H1
and H2.  H1 and H2 must be of the same dimensions unless one is a
scalar.  The result is a logical array of the same dimensions, where
each element is an element-wise equality result.

If one of H1 or H2 is scalar, scalar expansion is performed and the
result will match the dimensions of the array that is not scalar.

TF = NE(H1, H2) stores the result in a logical array of the same
dimensions.

See also HANDLE, HANDLE/EQ, HANDLE/GE, HANDLE/GT, HANDLE/LE, HANDLE/LT

Documentation for handle/ne
doc handle.ne
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : A, B
- OutputNames : TF
- DefiningClass : handle

####lt#### : _< (LT)   Less than relation for handles._
Description:
H1 < H2 performs element-wise comparisons between handle arrays H1 and
H2.  H1 and H2 must be of the same dimensions unless one is a scalar.
The result is a logical array of the same dimensions, where each
element is an element-wise < result.

If one of H1 or H2 is scalar, scalar expansion is performed and the
result will match the dimensions of the array that is not scalar.

TF = LT(H1, H2) stores the result in a logical array of the same
dimensions.

See also HANDLE, HANDLE/EQ, HANDLE/GE, HANDLE/GT, HANDLE/LE, HANDLE/NE

Documentation for handle/lt
doc handle.lt
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : A, B
- OutputNames : TF
- DefiningClass : handle

####gt#### : _> (GT)   Greater than relation for handles._
Description:
H1 > H2 performs element-wise comparisons between handle arrays H1 and
H2.  H1 and H2 must be of the same dimensions unless one is a scalar.
The result is a logical array of the same dimensions, where each
element is an element-wise > result.

If one of H1 or H2 is scalar, scalar expansion is performed and the
result will match the dimensions of the array that is not scalar.

TF = GT(H1, H2) stores the result in a logical array of the same
dimensions.

See also HANDLE, HANDLE/EQ, HANDLE/GE, HANDLE/LE, HANDLE/LT, HANDLE/NE

Documentation for handle/gt
doc handle.gt
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : A, B
- OutputNames : TF
- DefiningClass : handle

####le#### : _<= (LE)   Less than or equal relation for handles._
Description:
Handles are equal if they are handles for the same object.  All
comparisons use a number associated with each handle object.  Nothing
can be assumed about the result of a handle comparison except that the
repeated comparison of two handles in the same MATLAB session will
yield the same result.  The order of handle values is purely arbitrary
and has no connection to the state of the handle objects being
compared.

H1 <= H2 performs element-wise comparisons between handle arrays H1 and
H2.  H1 and H2 must be of the same dimensions unless one is a scalar.
The result is a logical array of the same dimensions, where each
element is an element-wise >= result.

If one of H1 or H2 is scalar, scalar expansion is performed and the
result will match the dimensions of the array that is not scalar.

TF = LE(H1, H2) stores the result in a logical array of the same
dimensions.

See also HANDLE, HANDLE/EQ, HANDLE/GE, HANDLE/GT, HANDLE/LT, HANDLE/NE

Documentation for handle/le
doc handle.le
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : A, B
- OutputNames : TF
- DefiningClass : handle

####ge#### : _>= (GE)   Greater than or equal relation for handles._
Description:
H1 >= H2 performs element-wise comparisons between handle arrays H1 and
H2.  H1 and H2 must be of the same dimensions unless one is a scalar.
The result is a logical array of the same dimensions, where each
element is an element-wise >= result.

If one of H1 or H2 is scalar, scalar expansion is performed and the
result will match the dimensions of the array that is not scalar.

TF = GE(H1, H2) stores the result in a logical array of the same
dimensions.

See also HANDLE, HANDLE/EQ, HANDLE/GT, HANDLE/LE, HANDLE/LT, HANDLE/NE

Documentation for handle/ge
doc handle.ge
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : A, B
- OutputNames : TF
- DefiningClass : handle

####delete#### : _Delete a handle object._
Description:
The DELETE method deletes a handle object but does not clear the handle
from the workspace.  A deleted handle is no longer valid.

DELETE(H) deletes the handle object H, where H is a scalar handle.

See also HANDLE, HANDLE/ISVALID, CLEAR

Documentation for handle/delete
doc handle.delete
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : obj
- OutputNames : [N/A]
- DefiningClass : handle

####isvalid#### : _Test handle validity._
Description:
TF = ISVALID(H) performs an element-wise check for validity on the
handle elements of H.  The result is a logical array of the same
dimensions as H, where each element is the element-wise validity
result.

A handle is invalid if it has been deleted or if it is an element
of a handle array and has not yet been initialized.

See also HANDLE, HANDLE/DELETE

Documentation for handle/isvalid
doc handle.isvalid
- Access : public
- Static : false
- Abstract : false
- Sealed : true
- ExplicitConversion : false
- Hidden : false
- InputNames : obj
- OutputNames : validity
- DefiningClass : handle

####findprop#### : _Find property of MATLAB handle object._
Description:
p = FINDPROP(H,PROPNAME) finds and returns the META.PROPERTY object
associated with property name PROPNAME of scalar handle object H.
PROPNAME can be a string scalar or character vector.  It can be the
name of a property defined by the class of H or a dynamic property
added to scalar object H.

If no property named PROPNAME exists for object H, an empty
META.PROPERTY array is returned.

See also HANDLE, HANDLE/FINDOBJ, DYNAMICPROPS, META.PROPERTY

Documentation for handle/findprop
doc handle.findprop
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : object, propname
- OutputNames : prop
- DefiningClass : handle

####notify#### : _Notify listeners of event._
Description:
NOTIFY(H, eventname) notifies listeners added to the event named
eventname for handle object array H that the event is taking place.
eventname can be a string scalar or character vector.
H is the array of handles to the event source objects, and 'eventname'
must be a character vector.

NOTIFY(H,eventname,ed) provides a way of encapsulating information
about an event which can then be accessed by each registered listener.
ed must belong to the EVENT.EVENTDATA class.

See also HANDLE, HANDLE/ADDLISTENER, HANDLE/LISTENER, EVENT.EVENTDATA, EVENTS

Documentation for handle/notify
doc handle.notify
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, eventname
- OutputNames : [N/A]
- DefiningClass : handle

####notify#### : _Notify listeners of event._
Description:
NOTIFY(H, eventname) notifies listeners added to the event named
eventname for handle object array H that the event is taking place.
eventname can be a string scalar or character vector.
H is the array of handles to the event source objects, and 'eventname'
must be a character vector.

NOTIFY(H,eventname,ed) provides a way of encapsulating information
about an event which can then be accessed by each registered listener.
ed must belong to the EVENT.EVENTDATA class.

See also HANDLE, HANDLE/ADDLISTENER, HANDLE/LISTENER, EVENT.EVENTDATA, EVENTS

Documentation for handle/notify
doc handle.notify
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, eventname, eventdata
- OutputNames : [N/A]
- DefiningClass : handle

####addlistener#### : _Add listener for event._
Description:
el = ADDLISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = ADDLISTENER(hSource, PropName, Eventname, Callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be
a string scalar or character vector.  PropName must be a single
property name specified as string scalar or character vector, or a
collection of property names specified as a cell array of character
vectors or a string array, or as an array of one or more
meta.property objects.  The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, addlistener returns an event.listener.  To remove a
listener, delete the object returned by addlistener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.

ADDLISTENER binds the listener's lifecycle to the object that is the
source of the event.  Unless you explicitly delete the listener, it is
destroyed only when the source object is destroyed.  To control the
lifecycle of the listener independently from the event source object,
use listener or the event.listener constructor to create the listener.

See also LISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/addlistener
doc handle.addlistener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, eventname, callback
- OutputNames : L
- DefiningClass : handle

####listener#### : _Add listener for event without binding the listener to the source object._
Description:
el = LISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = LISTENER(hSource, PropName, Eventname, callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be a
string sclar or character vector.  PropName must be either a single
property name specified as a string scalar or character vector, or
a collection of property names specified as a cell array of character
vectors or a string array, or as an array of one ore more
meta.property objects. The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, listener returns an event.listener.  To remove a
listener, delete the object returned by listener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.  Calling delete(el) on the listener
object deletes the listener, which means the event no longer causes
the callback function to execute.

LISTENER does not bind the listener's lifecycle to the object that is
the source of the event.  Destroying the source object does not impact
the lifecycle of the listener object.  A listener created with LISTENER
must be destroyed independently of the source object.  Calling
delete(el) explicitly destroys the listener. Redefining or clearing
the variable containing the listener can delete the listener if no
other references to it exist.  To tie the lifecycle of the listener to
the lifecycle of the source object, use addlistener.

See also ADDLISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/listener
doc handle.listener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, eventname, callback
- OutputNames : L
- DefiningClass : handle

####addlistener#### : _Add listener for event._
Description:
el = ADDLISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = ADDLISTENER(hSource, PropName, Eventname, Callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be
a string scalar or character vector.  PropName must be a single
property name specified as string scalar or character vector, or a
collection of property names specified as a cell array of character
vectors or a string array, or as an array of one or more
meta.property objects.  The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, addlistener returns an event.listener.  To remove a
listener, delete the object returned by addlistener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.

ADDLISTENER binds the listener's lifecycle to the object that is the
source of the event.  Unless you explicitly delete the listener, it is
destroyed only when the source object is destroyed.  To control the
lifecycle of the listener independently from the event source object,
use listener or the event.listener constructor to create the listener.

See also LISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/addlistener
doc handle.addlistener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, properties, eventname, callback
- OutputNames : L
- DefiningClass : handle

####listener#### : _Add listener for event without binding the listener to the source object._
Description:
el = LISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = LISTENER(hSource, PropName, Eventname, callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be a
string sclar or character vector.  PropName must be either a single
property name specified as a string scalar or character vector, or
a collection of property names specified as a cell array of character
vectors or a string array, or as an array of one ore more
meta.property objects. The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, listener returns an event.listener.  To remove a
listener, delete the object returned by listener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.  Calling delete(el) on the listener
object deletes the listener, which means the event no longer causes
the callback function to execute.

LISTENER does not bind the listener's lifecycle to the object that is
the source of the event.  Destroying the source object does not impact
the lifecycle of the listener object.  A listener created with LISTENER
must be destroyed independently of the source object.  Calling
delete(el) explicitly destroys the listener. Redefining or clearing
the variable containing the listener can delete the listener if no
other references to it exist.  To tie the lifecycle of the listener to
the lifecycle of the source object, use addlistener.

See also ADDLISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/listener
doc handle.listener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, properties, eventname, callback
- OutputNames : L
- DefiningClass : handle

####addlistener#### : _Add listener for event._
Description:
el = ADDLISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = ADDLISTENER(hSource, PropName, Eventname, Callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be
a string scalar or character vector.  PropName must be a single
property name specified as string scalar or character vector, or a
collection of property names specified as a cell array of character
vectors or a string array, or as an array of one or more
meta.property objects.  The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, addlistener returns an event.listener.  To remove a
listener, delete the object returned by addlistener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.

ADDLISTENER binds the listener's lifecycle to the object that is the
source of the event.  Unless you explicitly delete the listener, it is
destroyed only when the source object is destroyed.  To control the
lifecycle of the listener independently from the event source object,
use listener or the event.listener constructor to create the listener.

See also LISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/addlistener
doc handle.addlistener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, propertyname, eventname, callback
- OutputNames : L
- DefiningClass : handle

####listener#### : _Add listener for event without binding the listener to the source object._
Description:
el = LISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = LISTENER(hSource, PropName, Eventname, callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be a
string sclar or character vector.  PropName must be either a single
property name specified as a string scalar or character vector, or
a collection of property names specified as a cell array of character
vectors or a string array, or as an array of one ore more
meta.property objects. The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, listener returns an event.listener.  To remove a
listener, delete the object returned by listener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.  Calling delete(el) on the listener
object deletes the listener, which means the event no longer causes
the callback function to execute.

LISTENER does not bind the listener's lifecycle to the object that is
the source of the event.  Destroying the source object does not impact
the lifecycle of the listener object.  A listener created with LISTENER
must be destroyed independently of the source object.  Calling
delete(el) explicitly destroys the listener. Redefining or clearing
the variable containing the listener can delete the listener if no
other references to it exist.  To tie the lifecycle of the listener to
the lifecycle of the source object, use addlistener.

See also ADDLISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/listener
doc handle.listener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, propertyname, eventname, callback
- OutputNames : L
- DefiningClass : handle

####addlistener#### : _Add listener for event._
Description:
el = ADDLISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = ADDLISTENER(hSource, PropName, Eventname, Callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be
a string scalar or character vector.  PropName must be a single
property name specified as string scalar or character vector, or a
collection of property names specified as a cell array of character
vectors or a string array, or as an array of one or more
meta.property objects.  The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, addlistener returns an event.listener.  To remove a
listener, delete the object returned by addlistener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.

ADDLISTENER binds the listener's lifecycle to the object that is the
source of the event.  Unless you explicitly delete the listener, it is
destroyed only when the source object is destroyed.  To control the
lifecycle of the listener independently from the event source object,
use listener or the event.listener constructor to create the listener.

See also LISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/addlistener
doc handle.addlistener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, propertyname, eventname, callback
- OutputNames : L
- DefiningClass : handle

####listener#### : _Add listener for event without binding the listener to the source object._
Description:
el = LISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = LISTENER(hSource, PropName, Eventname, callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be a
string sclar or character vector.  PropName must be either a single
property name specified as a string scalar or character vector, or
a collection of property names specified as a cell array of character
vectors or a string array, or as an array of one ore more
meta.property objects. The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, listener returns an event.listener.  To remove a
listener, delete the object returned by listener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.  Calling delete(el) on the listener
object deletes the listener, which means the event no longer causes
the callback function to execute.

LISTENER does not bind the listener's lifecycle to the object that is
the source of the event.  Destroying the source object does not impact
the lifecycle of the listener object.  A listener created with LISTENER
must be destroyed independently of the source object.  Calling
delete(el) explicitly destroys the listener. Redefining or clearing
the variable containing the listener can delete the listener if no
other references to it exist.  To tie the lifecycle of the listener to
the lifecycle of the source object, use addlistener.

See also ADDLISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/listener
doc handle.listener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, propertyname, eventname, callback
- OutputNames : L
- DefiningClass : handle

####addlistener#### : _Add listener for event._
Description:
el = ADDLISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = ADDLISTENER(hSource, PropName, Eventname, Callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be
a string scalar or character vector.  PropName must be a single
property name specified as string scalar or character vector, or a
collection of property names specified as a cell array of character
vectors or a string array, or as an array of one or more
meta.property objects.  The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, addlistener returns an event.listener.  To remove a
listener, delete the object returned by addlistener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.

ADDLISTENER binds the listener's lifecycle to the object that is the
source of the event.  Unless you explicitly delete the listener, it is
destroyed only when the source object is destroyed.  To control the
lifecycle of the listener independently from the event source object,
use listener or the event.listener constructor to create the listener.

See also LISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/addlistener
doc handle.addlistener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, propertynames, eventname, callback
- OutputNames : L
- DefiningClass : handle

####listener#### : _Add listener for event without binding the listener to the source object._
Description:
el = LISTENER(hSource, Eventname, callbackFcn) creates a listener
for the event named Eventname.  The source of the event is the handle
object hSource.  If hSource is an array of source handles, the listener
responds to the named event on any handle in the array.  callbackFcn
is a function handle that is invoked when the event is triggered.

el = LISTENER(hSource, PropName, Eventname, callback) adds a
listener for a property event.  Eventname must be one of
'PreGet', 'PostGet', 'PreSet', or 'PostSet'. Eventname can be a
string sclar or character vector.  PropName must be either a single
property name specified as a string scalar or character vector, or
a collection of property names specified as a cell array of character
vectors or a string array, or as an array of one ore more
meta.property objects. The properties must belong to the class of
hSource.  If hSource is scalar, PropName can include dynamic
properties.

For all forms, listener returns an event.listener.  To remove a
listener, delete the object returned by listener.  For example,
delete(el) calls the handle class delete method to remove the listener
and delete it from the workspace.  Calling delete(el) on the listener
object deletes the listener, which means the event no longer causes
the callback function to execute.

LISTENER does not bind the listener's lifecycle to the object that is
the source of the event.  Destroying the source object does not impact
the lifecycle of the listener object.  A listener created with LISTENER
must be destroyed independently of the source object.  Calling
delete(el) explicitly destroys the listener. Redefining or clearing
the variable containing the listener can delete the listener if no
other references to it exist.  To tie the lifecycle of the listener to
the lifecycle of the source object, use addlistener.

See also ADDLISTENER, EVENT.LISTENER, HANDLE, NOTIFY, DELETE, META.PROPERTY, EVENTS

Documentation for handle/listener
doc handle.listener
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : sources, propertynames, eventname, callback
- OutputNames : L
- DefiningClass : handle

####findobj#### : _Find objects matching specified conditions._
Description:
The FINDOBJ method of the HANDLE class follows the same syntax as the
MATLAB FINDOBJ command, except that the first argument must be an array
of handles to objects.

HM = FINDOBJ(H, <conditions>) searches the handle object array H and
returns an array of handle objects matching the specified conditions.
Only the public members of the objects of H are considered when
evaluating the conditions.

See also FINDOBJ, HANDLE

Documentation for handle/findobj
doc handle.findobj
- Access : public
- Static : false
- Abstract : false
- Sealed : false
- ExplicitConversion : false
- Hidden : false
- InputNames : H, varargin
- OutputNames : HM
- DefiningClass : handle
