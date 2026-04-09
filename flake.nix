{
  description = "Fortran environment with gnuplot";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
    in
    {
      packages.${system}.default = 
        pkgs.stdenv.mkDerivation {
          name = "programa";
          src = ./.;
          buildInputs = [ pkgs.gfortran ];
          buildPhase = ''
            gfortran -o programa programa.f subrutina.f -lm
          '';
          installPhase = ''
            mkdir -p $out
            cp programa cuadratica.dat $out/
          '';
        };

      devShells.${system}.default = pkgs.mkShell {
        buildInputs = with pkgs; [
          gfortran
          gnuplot
        ];
      };
    };
}
