# mustBeFile

**Filetype:** _MATLAB&reg; function_

**Synopsis:** _Validate value is a filename (or set of filenames) on the path_

MUSTBEFILE is a validation function
which wraps to isfile() and issues an error if the input argument is not comprised solely of
existing files on the path. (The input can be a string array, character array, or
cell array of character vectors (aka "cell-string") of any size.)

### Usage ###

[] = MUSTBEFILE( A )

### References ###

See also

<https://www.mathworks.com/help/matlab/ref/isfile.html isfile>

<https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html validation functions>


### Attributes


- nInputs : 1

- nOutputs : 0
