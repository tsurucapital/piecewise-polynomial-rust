{ pkgs ? import ./nix/nixpkgs.nix }:

with pkgs;

let
  rustChannel = rustChannelOf
    { rustToolchain = ./rust-toolchain;
      sha256 = "sha256:0qp2gc5dm92wzrh6b2mqajdd1lplpl16l5f7km7d6hyx82brm3ij";
      date = "2020-07-17";
    };
in mkShell {
  buildInputs = [
    rustChannel.rust
    # Dev tools
    cargo-asm linuxPackages.perf
  ];
  # Required by test-suite and in general let's set a uniform one.
  LANG="C.UTF-8";
}
