{
  description = "bench6";
  nixConfig.bash-prompt = "\[nix-develop\]$ ";

  inputs.jungle.url = "git+https://git.sr.ht/~rodarima/jungle";

  outputs = { self, jungle, ... }:
  let
    targetMachine = jungle.outputs.nixosConfigurations.hut;
    pkgs = targetMachine.pkgs;
  in {
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
          nanos6
          nodes
          nosv
          mpi
          tampi
        ];

        enableParallelBuilding = false;
        hardeningDisable = [ "all" ];
        dontStrip = true;
      };
    };
  };
}
