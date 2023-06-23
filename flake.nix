{
  description = "bench6";
  nixConfig.bash-prompt = "\[nix-develop\]$ ";

  inputs.jungle.url = "jungle";

  outputs = { self, jungle }:
  let pkgs = jungle.packages.x86_64-linux.hut; in
  {
    packages.x86_64-linux.default = pkgs.stdenv.mkDerivation rec {
      pname = "bench6";
      version = "local";

      src = self.outPath;

      buildInputs = with pkgs; [
        cmake
        bsc.clangOmpss2
        bsc.nanos6
        bsc.nodes
        bsc.mpi
        bsc.tampi
      ];

      enableParallelBuilding = false;
      cmakeFlags = [
        "-DCMAKE_C_COMPILER=clang"
        "-DCMAKE_CXX_COMPILER=clang++"
      ];
      hardeningDisable = [ "all" ];
      dontStrip = true;
    };
  };
}
