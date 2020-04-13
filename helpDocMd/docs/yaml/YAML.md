# YAML

**Filetype:** _MATLAB&reg; classdef_

**Synopsis:** _Serialize a matlab variable to yaml format_

[ X ] = YAML.load( S )
[ S ] = YAML.dump( X )

[ X ] = YAML.read( filepath )
YAML.write( filepath, X )

YAML.LOAD takes YAML string S and returns matlab variable X.
YAML.DUMP takes matlab variable X and converts to YAML string S.
YAML.READ and YAML.WRITE are convenient methods to load and dump
YAML format directly from a file.

Examples:
To serialize matlab object

    >> X = struct('matrix', rand(3,4), 'char', 'hello');
    >> S = YAML.dump(X);
    >> disp(S);
    matrix:
    - [0.9571669482429456, 0.14188633862721534]
    - [0.4853756487228412, 0.421761282626275]
    - [0.8002804688888001, 0.9157355251890671]
    char: hello

To decode yaml string

    >> X = YAML.load(S);
    >> disp(X)
      matrix: [3x2 double]
        char: 'hello'

See also: xmlread xmlwrite

    Documentation for YAML
       doc YAML


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


### JARFILE

**Synopsis:** _YAML.JARFILE is a property._

<table>
<table border=1><tr><th>Dependent</th><th>Constant</th><th>Abstract</th><th>Transient</th><th>Hidden</th><th>GetObservable</th><th>SetObservable</th><th>AbortSet</th><th>NonCopyable</th><th>HasDefault</th></tr>
<tr><td>false</td><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td><td>true</td><td>true</td></tr>
</table>

- GetAccess : public
- SetAccess : none
- PartialMatchPriority : [N/A] 
- GetMethod : 
- SetMethod : 
- DefaultValue : /Users/ryan/Projects/General/scripts/shim/helpDocMd/src/yaml/java/snakeyaml-1.9.jar
- Validation : [N/A] 
- DefiningClass : YAML

---
## Methods


---


### write

**Synopsis**: _serialize and write yaml data to file_ 

 WRITE serialize and write yaml data to file


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : filepath, X
- OutputNames : [N/A] 
- DefiningClass : YAML

---


### read

**Synopsis**: _read and decode yaml data from file_ 

 READ read and decode yaml data from file


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : filepath
- OutputNames : X
- DefiningClass : YAML

---


### dump

**Synopsis**: _serialize matlab object into yaml string_ 

 DUMP serialize matlab object into yaml string


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : X
- OutputNames : S
- DefiningClass : YAML

---


### load

**Synopsis**: _load matlab object from yaml string_ 

 LOAD load matlab object from yaml string


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : S
- OutputNames : X
- DefiningClass : YAML

---


### jarfile

**Synopsis**: _path to the SnakeYAML jar file_ 

 JARFILE path to the SnakeYAML jar file


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : [N/A] 
- OutputNames : S
- DefiningClass : YAML

---


### dump_data

**Synopsis**: _convert_ 

 DUMP_DATA convert


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : private
- InputNames : r
- OutputNames : result
- DefiningClass : YAML

---


### merge_cell

**Synopsis**: _convert cell array to native matrix_ 

 MERGE_CELL convert cell array to native matrix


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : private
- InputNames : r
- OutputNames : result
- DefiningClass : YAML

---


### load_data

**Synopsis**: _recursively convert java objects_ 

 LOAD_DATA recursively convert java objects


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>true</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : private
- InputNames : r
- OutputNames : result
- DefiningClass : YAML

---


### YAML

**Synopsis**: _Serialize a matlab variable to yaml format_ 

[ X ] = YAML.load( S )
[ S ] = YAML.dump( X )

[ X ] = YAML.read( filepath )
YAML.write( filepath, X )

YAML.LOAD takes YAML string S and returns matlab variable X.
YAML.DUMP takes matlab variable X and converts to YAML string S.
YAML.READ and YAML.WRITE are convenient methods to load and dump
YAML format directly from a file.

Examples:
To serialize matlab object

    >> X = struct('matrix', rand(3,4), 'char', 'hello');
    >> S = YAML.dump(X);
    >> disp(S);
    matrix:
    - [0.9571669482429456, 0.14188633862721534]
    - [0.4853756487228412, 0.421761282626275]
    - [0.8002804688888001, 0.9157355251890671]
    char: hello

To decode yaml string

    >> X = YAML.load(S);
    >> disp(X)
      matrix: [3x2 double]
        char: 'hello'

See also: xmlread xmlwrite


#### Attributes:

<table>
<table border=1><tr><th>Static</th><th>Abstract</th><th>Sealed</th><th>ExplicitConversion</th><th>Hidden</th></tr>
<tr><td>false</td><td>false</td><td>false</td><td>false</td><td>false</td></tr>
</table>

- Access : public
- InputNames : [N/A] 
- OutputNames : obj
- DefiningClass : YAML

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
- DefiningClass : YAML
