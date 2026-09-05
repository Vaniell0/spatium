// The classic spinning ASCII donut (a1k0n lineage), built on
// spatium::render::terminal_canvas.hpp -- a torus swept by (theta, phi),
// rotated by two angles (A, B), projected with a per-pixel z-buffer, lit
// by a fixed light direction baked into the L formula below.
//
// The original donut.c hardcodes an extra "/2" in its y-projection to
// compensate for terminal cells being taller than wide -- that hack is
// exactly the "square pixel" problem terminal_canvas.hpp's
// query_terminal_size() measures instead of guesses; see cell_aspect
// below.
#include <spatium/render/terminal_canvas.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <thread>
#include <vector>

using namespace spatium::render;

int main() {
    const double R1 = 1.0;                  // tube radius
    const double R2 = 2.0;                  // donut radius (center to tube center)
    const double K2 = 12.5;                 // viewer distance -- pulled back ~2.5x; K1
                                             // below scales with K2 so apparent size is
                                             // unchanged, only perspective distortion drops
    const double two_pi = 2.0 * std::numbers::pi;

    // Custom palette -- write_terminal_frame_ascii takes any ordered
    // dim->bright character ramp, this isn't tied to the library default.
    constexpr std::string_view palette = " .:-=+*#%@";

    TerminalCursorGuard cursor_guard;

    // Re-queried every frame below rather than once up front: a resize
    // while the demo is running (dragging the terminal window, tmux pane
    // split/resize) would otherwise leave it rendering at a stale size --
    // ioctl() is cheap enough at 20fps that polling it costs nothing worth
    // avoiding with a SIGWINCH handler.
    int width = 0, height = 0;
    std::vector<std::vector<double>> intensity;
    std::vector<double> zbuffer;
    std::string frame_buf;

    double A = 0.0, B = 0.0;
    while (true) {
        TerminalSize term = query_terminal_size({.cols = 80, .rows = 45, .cell_aspect = 0.5});
        int new_width = term.cols;
        int new_height = std::max(1, term.rows - 1); // leave the last line for the shell prompt
        if (new_width != width || new_height != height) {
            width = new_width;
            height = new_height;
            intensity.assign(height, std::vector<double>(width, -1.0));
            zbuffer.assign(static_cast<std::size_t>(width) * height, 0.0);
        }
        const double cell_aspect = term.cell_aspect;

        const double K1 = width * K2 * 3.0 / (8.0 * (R1 + R2)); // projection scale

        // The reference donut.c's 0.07/0.02 spacing was tuned for an
        // ~80-column terminal. On a wider one the same angular step covers
        // proportionally more on-screen pixels (K1 grows with width), so
        // consecutive samples land more than a cell apart and leave gaps a
        // viewer reads as "transparent" patches (the far wall of the tube
        // showing through). Scaling both spacings down by 80/width keeps
        // the on-screen distance between samples roughly constant
        // regardless of terminal size.
        const double theta_spacing = 0.07 * 80.0 / width;
        const double phi_spacing = 0.02 * 80.0 / width;

        std::fill(zbuffer.begin(), zbuffer.end(), 0.0);
        for (auto& row : intensity) std::fill(row.begin(), row.end(), -1.0);

        double cosA = std::cos(A), sinA = std::sin(A);
        double cosB = std::cos(B), sinB = std::sin(B);

        for (double theta = 0.0; theta < two_pi; theta += theta_spacing) {
            double costheta = std::cos(theta), sintheta = std::sin(theta);

            for (double phi = 0.0; phi < two_pi; phi += phi_spacing) {
                double cosphi = std::cos(phi), sinphi = std::sin(phi);

                double circlex = R2 + R1 * costheta;
                double circley = R1 * sintheta;

                double x = circlex * (cosB * cosphi + sinA * sinB * sinphi) - circley * cosA * sinB;
                double y = circlex * (sinB * cosphi - sinA * cosB * sinphi) + circley * cosA * cosB;
                double z = K2 + cosA * circlex * sinphi + circley * sinA;
                double ooz = 1.0 / z; // "one over z" -- also doubles as depth for the z-buffer

                int xp = static_cast<int>(width / 2.0 + K1 * ooz * x);
                // Measured cell_aspect replaces donut.c's hardcoded "/2".
                int yp = static_cast<int>(height / 2.0 - K1 * ooz * y * cell_aspect);

                // Surface-normal . light-direction, algebraically expanded
                // rather than formed as two explicit vectors. At A=B=0 this
                // reduces to sin(theta) - cos(theta)*sin(phi), i.e. the
                // torus's own (cos(theta)cos(phi), cos(theta)sin(phi),
                // sin(theta)) normal dotted with a FIXED light direction of
                // (0,-1,1) in the torus's own frame. The light doesn't move
                // -- the object rotates under it as A/B animate, which is
                // why the highlight sweeps across the surface as it spins.
                double L = cosphi * costheta * sinB - cosA * costheta * sinphi - sinA * sintheta +
                           cosB * (cosA * sintheta - costheta * sinA * sinphi);

                if (L <= 0.0) continue; // facing away from the light
                if (xp < 0 || xp >= width || yp < 0 || yp >= height) continue;

                std::size_t idx = static_cast<std::size_t>(yp) * width + xp;
                if (ooz > zbuffer[idx]) {
                    zbuffer[idx] = ooz;
                    intensity[yp][xp] = L / std::numbers::sqrt2; // L's range is [-sqrt2, sqrt2]
                }
            }
        }

        frame_buf.clear();
        write_terminal_frame_ascii(frame_buf, intensity, palette);
        std::fwrite(frame_buf.data(), 1, frame_buf.size(), stdout);
        std::fflush(stdout);

        A += 0.02; // half the rotation speed, plus a longer frame delay below
        B += 0.01;
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
}
