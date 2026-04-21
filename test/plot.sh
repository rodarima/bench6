#!/bin/sh

#repo='https://gitlab.pm.bsc.es/rarias/bench6#plot' # Remote
repo='.#plot' # Local

nix develop "$repo" -c sh -c 'Rscript $self/plot/miniapps.R test/data/miniapps.csv'
