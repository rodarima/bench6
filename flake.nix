{
  description = "bench6";
  nixConfig.bash-prompt = "\[nix-develop\]$ ";

  inputs.jungle.url = "path:/home/Computational/rarias/jungle";
  inputs.jungle.inputs.bscpkgs.url = "path:/home/Computational/rarias/bscpkgs";

  outputs = { self, jungle, ... }:
  let
    targetMachine = jungle.outputs.nixosConfigurations.hut;
    pkgs = targetMachine.pkgs;
  in {
    packages.x86_64-linux.default = pkgs.stdenv.mkDerivation rec {
      pname = "bench6";
      version = if self ? shortRev then self.shortRev else "dirty";

      src = self.outPath;

      buildInputs = with pkgs; [
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
}
