# mustBeStringOrCharOrCellstr

**Filetype:** _MATLAB&reg; function_

**Synopsis:** _Validate value is string, char, or cell of chars_

[] = mustBeStringOrCharOrCellstr( A )

MUSTBESTRINGORCHARORCELLSTR is a validation function which issues an error if
the input argument is not a string or character array, i.e.:
if ~( isstring( A ) || ischar( A ) || iscellstr( A ) )


### Attributes


- nInputs : 1

- nOutputs : 0
