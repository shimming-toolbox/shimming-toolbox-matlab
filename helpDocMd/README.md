# helpDocMd

`functions` to `help` extract ***.m***bedded source `doc; +`export to ***.md***

- [Overview](#overview)
- [Getting started](#getting-started)
- [Basics](#basics)
- [References](#references)
- [FAQe](#faqe)
- [License](#license)

## Overview

**_helpDocMd_** aims to provide an easy pipeline for publishing custom
MATLAB&reg; documentation:

+ particularly class-based API documentation  

+ specifically as [Markdown] to be readily hosted online ([1][mkdocs], [2][github], [3][readthedocs]). 

[mkdocs]: https://www.mkdocs.org
[github]: https://pages.github.com/
[readthedocs]: https://readthedocs.org/
[Markdown]: https://daringfireball.net/projects/markdown/

## Getting started

### Dependencies

- MATLAB&reg; Version R2019b or later.

### Installation

1. Clone this repository or download the source code.

2. Start up Matlab and add the downloaded 'src' directory along with its
   subfolders to the path by typing `addpath(genpath( '.../src/' )  ;
   which('Documentor')` into the command prompt (replacing the ellipsis with
   the path to the downloaded folder). The filepath to Documentor.m should be
   displayed.

## Basics

### Printing

The first and most difficult step to publishing your source code documentation:
*Write it!* (specifically, the
[header comment](https://www.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html)).

To publish it, the process is simple:

1. Assemble the list of .m files and/or code-containing directories that you
   want to include as entries in a string vector: let's call it `src`.

2. Initialize a `Documentor` object (e.g. one called `Dr`) with your list:  

    Dr = Documentor( src ) ;  

3. To print all of the documentation to Markdown text files, call:  

    Dr.printdoc( ) ;  

#### Example

To print the documentation for the Documentor class itself:

    Dr = Documentor( which('Documentor') ) ;

    Dr.printdoc( ) ;

If successful, the documentation file path is displayed in the command window.
Once formated to HTML on GitHub, it looks like
[this](https://github.com/neuropoly/realtime_shimming/blob/helpDocMd/helpDocMd/docs/Documentor.md).

## FAQe 

—*A fake FAQ:*

>What's wrong with the built-in `help` and `doc` commands?

Both are very useful, and the recently introduced 
[Live Editor](https://www.mathworks.com/help/matlab/matlab_prog/what-is-a-live-script-or-function.html)
has great publishing features for scripts, but as of 2020(a) there still isn't
any built-in functionality to *export* class documentation.

>Why not use alternative software, like
>[doxygen](https://github.com/simgunz/doxymatlab) or
>[sphinx](https://pypi.org/project/sphinxcontrib-matlabdomain/)?

There are advantages—simplicity, +potentially
[automating](https://www.mathworks.com/help/matlab/ref/meta.class-class.html)
some of the documentation itself—to keeping the software and the commenting
[syntax](https://www.mathworks.com/matlabcentral/answers/help/markup/)
"[in-house](https://youtu.be/LYmvsFqu8kg?t=35)"; particularly if we can assume
that the built-in publishing features will eventually be extended to classes.
(Albeit, the `doc` command already generates the—yet private—HTML that users
have been requesting for years ([4][stack1], [5][stack2]) so who knows...)

[stack1]: https://stackoverflow.com/questions/26242145/what-is-the-mathworks-way-to-generate-matlab-html-documentation
[stack2]: https://stackoverflow.com/questions/37562403/publish-matlab-class-documentation-to-html?rq=1

## References

- [test a Markdown sample](https://daringfireball.net/projects/markdown/dingus) for how it will display once reformatted to HTML.

Re: hosting documentation online:
- [MkDocs](https://www.mkdocs.org/)
- [Github](https://pages.github.com/)
- [ReadTheDocs](https://docs.readthedocs.io/en/stable/)

## License

Same for parent repo: [realtime_shimming](https://github.com/neuropoly/realtime_shimming/blob/master/LICENSE)
