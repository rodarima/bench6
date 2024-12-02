#!/bin/sh

# Running in owl doesn't work yet
#repo="$CI_PROJECT_URL"
#branch="$CI_COMMIT_REF_NAME"
#commit="$CI_COMMIT_SHA"
#url="$repo?ref=$branch&rev=$commit"
#srun --exclusive nix develop "$url" -c '$bench6src/test/bench.sh'

nix develop -c test/bench.sh
