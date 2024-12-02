{
  description = "bench6";
  nixConfig.bash-prompt = "\[nix-develop\]$ ";

  inputs.bscpkgs.url = "git+https://git.sr.ht/~rodarima/bscpkgs";
  #inputs.jungle.url = "git+https://git.sr.ht/~rodarima/jungle";
  #inputs.jungle.inputs.bscpkgs.follow = "bscpkgs";

  outputs = { self, bscpkgs, ... }:
  let
    #targetMachine = jungle.outputs.nixosConfigurations.hut;
    #pkgs = targetMachine.pkgs;

    pkgs = import bscpkgs.inputs.nixpkgs {
      system = "x86_64-linux";
      overlays = [ bscpkgs.bscOverlay ];
    };
  in {
    devShells.x86_64-linux.default = self.outputs.devShells.x86_64-linux.run;
    devShells.x86_64-linux.run = pkgs.mkShell {
      buildInputs = [
        self.outputs.packages.x86_64-linux.bench6
        pkgs.bigotes
      ];
      bench6 = self.outputs.packages.x86_64-linux.bench6;
      bench6src = self.outputs.packages.x86_64-linux.bench6.src;
    };
    packages.x86_64-linux = rec {
      default = bench6;
      bench6 = pkgs.stdenv.mkDerivation rec {
        pname = "bench6";
        version = if self ? shortRev then self.shortRev else "dirty";

        src = self.outPath;

        buildInputs = with pkgs; [
          bigotes
          cmake
          clangOmpss2
          openmp
          openmpv
          nanos6
          nodes
          nosv
          #mpi
          #tampi
        ];

        buildFlags = [ "VERBOSE=1" ];
        enableParallelBuilding = false;
        hardeningDisable = [ "all" ];
        dontStrip = true;
      };

      bench6Master = bench6.overrideAttrs (old: {
        src = builtins.fetchGit {
          url = "https://pm.bsc.es/gitlab/rarias/bench6.git";
          #ref = "master";
          # FIXME: Just for testing
          ref = "heat";
          rev = "466e50e511e11087f6e9298b3ac851fcab7d459c";
        };
      });
    };
  };
}
