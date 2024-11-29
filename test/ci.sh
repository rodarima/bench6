#!/bin/sh

env

srun --exclusive nix develop -c test/bench.sh
