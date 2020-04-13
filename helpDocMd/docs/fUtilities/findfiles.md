# findfiles

**Filetype:** _MATLAB&reg; function_

**Synopsis:** _Search a directory for filenames matching a pattern_

      
     [paths, List] = findfiles( sFolder, sPattern, isRecursive, isExcludingHidden, returnType )  

Looks for files and/or subfolders in `sFolder` with names matching `sPattern`
by calling Matlab function [dir] and returns the file paths as elements of a
string column vector `paths`. The single-element structs output by `dir()`
are arrayed and returned as `List.

__INPUTS__
      
    sFolder=["."]  
      The base directory of the search.    

    sPattern=["*.*"]  
      The searchPattern of interest. If provided as a string array, patterns
      are searched successively. (The default corresponds to including  all
      files with explicit with explicit file extensions.)

    isRecursive=[true|1]  
      Toggle to include (1) or exclude (0) subdirectories in the search.

    isExcludingHidden=[true|1]  
      Toggle to include (1) or exclude (0) hidden files (i.e. for
      Unix: filenames beginning with ".")

    returnType=["files"]  
      Selects what type of path elements are retained in the two outputs:  
      Options are: "files", "folders", or "both".

ETC

    For more info, refer to the documentation for  
    [dir](https://www.mathworks.com/help/matlab/ref/dir.html)  

See also DIR

    Other functions named findfiles

       Documentor/findfiles


### Attributes


- nInputs : 5

- nOutputs : 2
