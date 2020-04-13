# mustBeFileOrFolder

**Filetype:** _MATLAB&reg; function_

**Synopsis:** _Assert input consists exclusively of valid file system paths_

     
     [] = mustBeFileorFolder( A ) 

Throws an error if any single element of `A` is neither an existing file nor
folder path: namely, if `all( [ isfile( A ) | isfolder( A ) ] )` evaluates to
`false`.

__INPUTS__

    `A` 
      A string- or cellstr-array of any size, or a single-row character vector.

__ETC__

-[Validation functions](https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html)
-[isfolder](https://www.mathworks.com/help/matlab/ref/isfolder.html)
-[isfile][https://www.mathworks.com/help/matlab/ref/isfile.html)

See also
ISFOLDER, ISFILE


### Attributes


- nInputs : 1

- nOutputs : 0
