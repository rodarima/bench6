#!/bin/sh

set -e
set -x

# Allow impure evaluation so we fetch the latest commit from the repo
bench6_ref=$(nix build --print-out-paths --impure ".#bench6Master")
bench6_cur=$(nix build --print-out-paths ".#bench6")

# Add bigotes to the path
bigotes=$(nix build --print-out-paths 'jungle#bigotes')
export PATH="$bigotes/bin:$PATH"

bigotes "${bench6_ref}/bin/b6_heat_nanos6" -s 2048 -t 10 -b 64
bigotes "${bench6_cur}/bin/b6_heat_nanos6" -s 2048 -t 10 -b 64
