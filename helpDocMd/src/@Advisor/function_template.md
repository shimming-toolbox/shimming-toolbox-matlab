FUNCTION_NAME *Brief summary goes here* (Note 1)

__SYNTAX__ <--*MATLAB convention for the 1st header title.* (Note 2)   
>1 empty line +4x spaces from margin formats the following `as code` in [Markdown]. (Note 3)

    [out1] = FUNCTION_NAME( mandatoryInput1, mandIn2_A )  
    [out1] = FUNCTION_NAME( mandatoryInput1, mandIn2_B )  
    [out1] = FUNCTION_NAME( mandatoryInput1, ..., optionalInput1 )  
    [out2] = FUNCTION_NAME( mandatoryInput1, ..., optionalInput1, optionalInput2 )  

__DESCRIPTION__ <--*(2x underscores format to bold in Markdown, but italics in MATLAB Markup...)*
>Describe the *basic* functionality (what will the anticipated user really need to know). You should organize the function call and the description such that the most common use scenario can be addressed first. More exotic applications can then be alluded to, here, and then described in more detail after going over all the basics (Note 4). Example, for example:


__INPUTS__  
>Be mindful of character limits/line. leave 2 spaces after the last word on a line to enforce a line break, otherwise the 2 lines will be merged when converted to Markdown.  

   `mandatoryInput1`  
     The file path to the .jpeg image to be processed (as a string scalar or character vector)

   `mandIn2_A`  
     The name of the spatial kernel to be applied to `mandatoryInput1`. Options are: "Gaussian" or "median"

   `mandIn_C`  
     The function handle to the custom filtering function

   `optionalInput1=[DEFAULT1]`  
     The kernel size as a 2-element double-vector in units of pixels

   `optionalInput2=[DEFAULT2]`  
     Toggles whether noise statistics are included in the output.

__OUTPUTS__

__NOTES__

1. a) MATLAB convention is to use ALL CAPS for every instance of the function
name (that is, despite the fact common naming conventions generally insist against naming functions with capital letters, or using lowercase exclusively.
If a user queries `help` while running the MATLAB GUI, FUNCTION_NAME is
automatically reformatted to its proper case should it be different from the file/function is actually called "function_Name.m". However, except for maybe the first (title) instance, I would recommend not
adopting this convention since Matlab will generate the appropriate hyperlinks either way when working in the IDE, but it stands to mislead anyone calling `help` while running MATLAB in the terminal since the name won't be reformatted.

  b) The leading "H1" line is special as it is searchable via the [lookfor] function, so
  make sure it contains relevant terms! Also, keep it brief (<80 char including the function name) or it risks being truncated when displayed by `lookfor`.

2. MATLAB-style Markup and Markdown are similar but feature a number of differences. Assuming either is relevant to your case, you will need to choose which is best. This template is designed as a sort of compromise to have something that looks ok when `help` is called via the commmand line, while nevertheless being easily converted to Markdown.

[H1](https://www.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html)
[lookfor](https://www.mathworks.com/help/matlab/ref/lookfor.html)

3.
[Markdown](https://daringfireball.net/projects/markdown/syntax)
[UNIX](https://en.wikipedia.org/wiki/Man_page) favors "SYNOPSIS", others use "USAGE"

4. A rule of thumb: The total number of lines of code should be roughly comparable to the lines of text dedicated to documenting a "Function" (that is, a "feature" of some form: not a single .m file necessarily, but some functional component of software designated for use by other people)

This is not to suggest that it be overly detailed: to the contrary, it should be as concise as possible or it risks not being read and, by extension, not being used at all. -- detailing only the essentials and/or likely points of interest for a given user.
Rather, if the number of lines of documentation vs. lines of code for a typical function is not on the same order, more than likely your code is too complicated and should be refactored!

__ETC__

For more info, refer to the documentation for

See also
LOOKFOR
