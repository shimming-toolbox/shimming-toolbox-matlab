# helpDocMd

Utilities for documenting custom Matlab&reg source code.

- [Purpose](#purpose)
- [Getting started](#gettingstarted)
- [Basics](#basics)
- [References](#references)
- [License](#license)
 
## Purpose statement ~of the Obvious~

The most convenient way to document code manually is, of course, to embed
descriptive comments within the source code itself. That said, MATLAB&reg; is,
after all, a programming _language_: to the extent the code works, it should
_literally_ describe itself. The internal documentation viewer (via `doc` and
`help` commands) is, of course, very useful; however, there are occasions when
the displayed content is curiously stingy. (E.g., if a function can only accept
values of a certain type, clearly, in most cases this would deserve mention. Of
course, you can write it in manually into the file's header comment (and,
again, whenever the code changes...) but, not only is that redundant (it's
already in the *code*), it betrays the basic reason for writing code in the
first place: namely, to assign *your* menial work to a machine.

The purpose of this codebase is to guide, structure, and (to the extent
possible) autogenerate Matlab&reg; source documentation that is readily hosted
online (e.g. [1][mkdocs],[2][github],[3][readthedocs]) in simple, readable [Markdown][markdown].
Further, this should not require a peculiar tagging syntax to MATLAB's own
style of [markup][markup] or additional dedependecies ([e.g.][sphinx]).


[mkdocs]: https://www.mkdocs.org
[github]: https://pages.github.com/
[readthedocs]: https://readthedocs.org/
[markup]: https://www.mathworks.com/matlabcentral/answers/help/markup/
[markdown]: https://daringfireball.net/projects/markdown/ 
[sphinx]: https://pypi.org/project/sphinxcontrib-matlabdomain/ 

 
## Getting started

### Dependencies
 
- MATLAB&reg; Version R2019b or later. 
9.7.0.1216025 () Update 1 PR2019b

### Installation

1. Clone this repository or download the source code.

2. Start up Matlab and add the '/src/' folder to the path recursively to include its subfolders:
`addpath(genpath( '.../src/' )  ;`

## Basics
 
See the documentation for the `Documentor` class for examples. WIP!...

## References
 
## License

