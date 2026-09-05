{
  description = "Spatium — C++23 mathematical spaces library with Vulkan visualization";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    imgui-src = {
      url = "github:ocornut/imgui/v1.91.4";
      flake = false;
    };
  };

  outputs = { self, nixpkgs, imgui-src }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems
        (system: f (import nixpkgs { inherit system; }));

      mkProject = pkgs:
        let
          stdenv  = pkgs.gcc15Stdenv;
          lib     = pkgs.lib;
          isLinux = pkgs.stdenv.isLinux;

          src = lib.cleanSourceWith {
            src = ./.;
            filter = path: type:
              let b = builtins.baseNameOf path; in
              !(b == "build" || b == "result" || b == ".git" || b == ".direnv");
          };

          nativeDeps = with pkgs; [ cmake ninja pkg-config ]
            ++ lib.optionals isLinux [ pkgs.makeWrapper ];

          coreDeps = with pkgs; [ catch2_3 boost gbenchmark eigen ];

          gfxDeps = lib.optionals isLinux (with pkgs; [
            glfw3 vulkan-loader vulkan-headers vulkan-validation-layers shaderc
          ]);

          spatium = stdenv.mkDerivation {
            pname   = "spatium";
            version = "1.0.0";
            inherit src;

            nativeBuildInputs = nativeDeps;
            buildInputs       = coreDeps ++ gfxDeps;

            cmakeFlags = [
              "-DSPATIUM_BUILD_TESTS=OFF"
              "-DSPATIUM_BUILD_BENCHMARKS=ON"
              "-DSPATIUM_EIGEN=ON"
              "-DCMAKE_SKIP_BUILD_RPATH=ON"
              "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON"
              "-DIMGUI_DIR=${imgui-src}"
            ];

            # installPhase через CMake install(TARGETS), fallback на cp
            # Each example binary becomes spatium-<name>.  Vulkan-only ones
            # may be skipped when SPATIUM_BUILD_VIEWER cannot find its deps.
            installPhase = ''
              runHook preInstall
              mkdir -p $out/bin $out/lib $out/include
              install_demo() {
                local src=$1
                local name=$2
                if [ -f "$src" ]; then
                  install -Dm755 "$src" "$out/bin/spatium-$name"
                fi
              }
              install_demo examples/geometry_demo              geometry
              install_demo examples/showcase                   showcase
              install_demo examples/sets_demo                  sets
              install_demo examples/parametric_analytical_demo parametric
              install_demo examples/blackhole_demo             blackhole
              install_demo examples/blackhole_gr_demo          blackhole-gr
              install_demo examples/geodesic_curvature_grid_demo curvature-grid
              install_demo examples/geodesic_procgen_demo       geodesic-procgen
              install_demo examples/hyperbolic_tessellation_demo hyperbolic-tessellation
              install_demo examples/wormhole_demo               wormhole
              install_demo examples/tumbling_body_demo          tumbling
              install_demo examples/wave_ca_demo               wave
              install_demo examples/collatz_demo               collatz
              install_demo examples/burning_ship_demo          burning-ship
              install_demo examples/atom_demo                  atom-demo
              install_demo examples/primitives_demo            primitives-demo
              install_demo examples/terminal_donut_demo        terminal-donut
              [ -f benchmarks/spatium_bench ] && install -Dm755 benchmarks/spatium_bench $out/bin/spatium-bench
              for so in lib*.so; do [ -f "$so" ] && install -Dm755 "$so" $out/lib/; done
              cp -r ${src}/include/spatium $out/include/
              runHook postInstall
            '' + lib.optionalString isLinux ''
              for prog in spatium-atom-demo spatium-primitives-demo; do
                [ -f $out/bin/$prog ] && wrapProgram $out/bin/$prog \
                  --prefix VK_LAYER_PATH : "${pkgs.vulkan-validation-layers}/share/vulkan/explicit_layer.d"
              done
            '';

            meta = {
              description = "C++23 mathematical spaces library";
              license     = lib.licenses.mit;
              mainProgram = "spatium-primitives-demo";
              platforms   = lib.platforms.unix;
            };
          };

          mkApp = bin: { type = "app"; program = "${spatium}/bin/${bin}"; };
        in { inherit spatium stdenv mkApp isLinux; };

    in {
      packages = forAllSystems (pkgs: {
        default = (mkProject pkgs).spatium;
      });

      apps = forAllSystems (pkgs:
        let p = mkProject pkgs; in {
          default       = p.mkApp "spatium-primitives-demo";
          atom          = p.mkApp "spatium-atom-demo";
          primitives    = p.mkApp "spatium-primitives-demo";
          showcase      = p.mkApp "spatium-showcase";
          sets          = p.mkApp "spatium-sets";
          geometry      = p.mkApp "spatium-geometry";
          parametric    = p.mkApp "spatium-parametric";
          blackhole     = p.mkApp "spatium-blackhole";
          blackhole-gr  = p.mkApp "spatium-blackhole-gr";
          curvature-grid = p.mkApp "spatium-curvature-grid";
          geodesic-procgen = p.mkApp "spatium-geodesic-procgen";
          hyperbolic-tessellation = p.mkApp "spatium-hyperbolic-tessellation";
          wormhole      = p.mkApp "spatium-wormhole";
          tumbling      = p.mkApp "spatium-tumbling";
          wave          = p.mkApp "spatium-wave";
          collatz       = p.mkApp "spatium-collatz";
          burning-ship  = p.mkApp "spatium-burning-ship";
          terminal-donut = p.mkApp "spatium-terminal-donut";
          bench         = p.mkApp "spatium-bench";
        });

      devShells = forAllSystems (pkgs:
        let
          p = mkProject pkgs;
          # cudaPackages is unfree; allowUnfree is scoped to this second
          # nixpkgs instantiation instead of the flake's main `pkgs`, so
          # every other output (default shell, packages, apps) is
          # unaffected.
          pkgsUnfree = import nixpkgs { inherit (pkgs) system; config.allowUnfree = true; };
        in {
          default = (pkgs.mkShell.override { inherit (p) stdenv; }) {
            inputsFrom = [ p.spatium ];
            packages   = with pkgs; [ ccache clang-tools gdb ];
            shellHook = ''
              export IMGUI_DIR="${imgui-src}"
              export CCACHE_DIR="''${CCACHE_DIR:-$HOME/.cache/ccache}"
              export CMAKE_C_COMPILER_LAUNCHER=ccache
              export CMAKE_CXX_COMPILER_LAUNCHER=ccache
            '' + nixpkgs.lib.optionalString p.isLinux ''
              export VK_LAYER_PATH="${pkgs.vulkan-validation-layers}/share/vulkan/explicit_layer.d"
            '';
          };

          # Opt-in, Linux-only: nvcc + the CUDA runtime headers/libs for
          # the gpu/ geodesic kernel (-DSPATIUM_CUDA=ON). Not part of the
          # default shell -- cudaPackages is unfree and this shell is
          # only useful on a machine with an actual NVIDIA GPU driver
          # (nvcc alone will cross-compile without one, but running
          # anything needs it). `nix develop .#cuda`.
          cuda = if p.isLinux then
            (pkgsUnfree.mkShell.override { inherit (p) stdenv; }) {
              inputsFrom = [ p.spatium ];
              packages = with pkgsUnfree; [
                ccache clang-tools gdb
                cudaPackages.cudatoolkit cudaPackages.cuda_nvcc
              ];
              shellHook = ''
                export IMGUI_DIR="${imgui-src}"
                export VK_LAYER_PATH="${pkgs.vulkan-validation-layers}/share/vulkan/explicit_layer.d"
              '';
            }
          else
            pkgs.mkShell {
              shellHook = ''
                echo "spatium's CUDA devShell is Linux-only (this is ${pkgs.system})."
                exit 1
              '';
            };
        });

      formatter = forAllSystems (pkgs: pkgs.nixpkgs-fmt);
    };
}
