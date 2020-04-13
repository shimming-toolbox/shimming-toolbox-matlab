# abspath

**Filetype:** _MATLAB&reg; function_

**Synopsis:** _Validate and convert file system paths_

     
     [pathOut, pathType] = abspath( pathIn )  

Wraps to the Matlab function `fileattrib` to check each element of the input
array `pathIn` for valid file system paths (which can be relative or
full/absolute). It returns two arrays of the same type (string, char, or
cellstr) and size (arbitrary) as the input.

The possibilities for a given return element `(i)` are:

`pathIn(i)`  | `pathOut(i)`| `pathType(i)`
-------------|-------------|-------------------------------------
is a file   | full path   | "file"
is a folder | full path   | "directory"
else        |   ""        | the error message from [fileattrib]

**Note** for invalid input paths, the error message returned from
`fileattrib( pathIn(i) )` will presumably be 'No such file or directory.'
However the MATLAB documentation is unclear whether other possibilities exist
for the case where the function is called with a single argument, as is done
here. For more info, refer to the MATLAB documentation:
[fleattrib](https://www.mathworks.com/help/matlab/ref/fileattrib.html/)

See also
FILEATTRIB


### Attributes


- nInputs : 1

- nOutputs : 2
