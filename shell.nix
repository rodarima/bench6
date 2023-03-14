let
  pkgs = import (builtins.fetchTarball
    "https://pm.bsc.es/gitlab/rarias/bscpkgs/-/archive/master/bscpkgs-master.tar.gz");

  rWrapper = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [ tidyverse rjson jsonlite egg viridis ];
  };
in
  pkgs.mkShell {
    nativeBuildInputs = [
      pkgs.bsc.clangOmpss2
      pkgs.bsc.nanos6
      rWrapper
    ];
  }
