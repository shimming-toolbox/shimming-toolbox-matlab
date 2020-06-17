#!/bin/sh

set -e

# currently in Shimming-toolbox repo named "s"
git branch
cd ..
ls
# Clone doc generation software
#git clone https://github.com/shimming-toolbox/helpDocMd.git

# Build API documentation, Matlab (.m) -> Markdown (.md)
#   Beware: you MUST put `quit` at the end; Matlab does not auto-quit.
#   https://www.mathworks.com/matlabcentral/answers/523194-matlab-script-in-batch-from-unix-command-line
# TODO: broken.
#  1. writes to temp/shimming-toolbox/docs instead of docs/
#  2. incompatible with Matlab R2019a which is what's on our CI machine.
/usr/local/MATLAB/R2020a/bin/matlab -nodisplay -nosplash -r "run('./s/generate_doc.m');exit"
cd s
# run mkdocs
# ported from https://github.com/mhausenblas/mkdocs-deploy-gh-pages/blob/master/action.sh
# with hints from https://github.com/DavidS/jekyll-deploy/blob/master/entrypoint.rb and https://docs.travis-ci.com/user/deployment/pages/
# maybe it should be its own script?

# install mkdocs
# TODO: move this to a dependencies phase instead?
#export PATH="~/.local/bin":"$PATH"
#pip install --user mkdocs

# we need to do this in an explicit venv since our build agent
# isn't fancy enough to spawn fresh containers/VMs for us
pip install --user virtualenv
VENV=$(mktemp -d)
trap 'rm -rf $VENV' EXIT  # cleanup after ourselves
python -m virtualenv "$VENV"
. "$VENV"/bin/activate
pip install --no-cache-dir --ignore-installed mkdocs

# Make sure GH_PAGES_TOKEN is set up in pipeline
if [ -n "${GH_PAGES_TOKEN}" ]; then
    # XXX this is destructive to the repo, but we're assuming this is being run in CI so no matter?

    #git config http."https://github.com/".extraheader "Authorization: token $GH_PAGES_TOKEN" # authenticate over HTTPS with the given token
    # ^ this doesn't work. github will accept 'token ...' for API calls, but not for git+https://
    # so instead we have to generate an HTTP basic auth:
    #auth="$(echo -n "x-access-token:$GH_PAGES_TOKEN" | base64 | tr -d '\r' | tr -d '\n')"
    #git config http."https://github.com/".extraheader "Authorization: basic "$auth"" # authenticate over HTTPS with the given token

    git remote set-url --push origin $(git remote get-url origin | sed 's!https://!https://x-access-token:'"$GH_PAGES_TOKEN"'@!')

    git config url."https://github".insteadOf "git@github.com:"  # force HTTPS, not SSH

    # if you can parse out the remote URL, you could instead try:
    # remote_repo="https://x-access-token:${GH_PAGES_TOKEN}@github.com/${GITHUB_REPOSITORY}.git" #
    # git remote rm origin
    # git remote add origin "${remote_repo}"
fi

if [ -n "${CUSTOM_DOMAIN}" ]; then
    echo "${CUSTOM_DOMAIN}" > "docs/CNAME"
fi


mkdocs gh-deploy --config-file "./mkdocs.yml" --force

