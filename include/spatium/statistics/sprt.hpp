#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::statistics {

// ── Whether one player is stronger than another, from as few games as the evidence needs ─
//
// Wald's sequential probability ratio test on games won and lost: H0 says the two are equal (a
// win is as likely as a loss), H1 that the first is stronger by `elo` points (a win has the chance
// 1 / (1 + 10^(-elo / 400)) against a loss). After each game the log-likelihood ratio
//
//   llr = wins log(p1 / p0) + losses log((1 - p1) / (1 - p0)),   p0 = 1/2
//
// is held against two bounds: above log((1 - beta) / alpha) H1 is taken, below log(beta / (1 -
// alpha)) H0, between them the games go on. So an equal player is taken as stronger with chance
// about alpha, a player `elo` stronger is missed with chance about beta, and the number of games
// is whatever it takes, few when the two are far apart. A draw says nothing about which is the
// stronger under this model and is left out.
//
// A fixed number of games scored against a threshold has no such guarantee: 24 games at a score
// of 0.55 take an equal player as stronger once in four.

struct Sprt {
    double elo = 35;
    double alpha = 0.05, beta = 0.05;

    // The chance of a win against a loss for a player `elo` points the stronger.
    static double win_chance(double elo) { return 1 / (1 + std::pow(10.0, -elo / 400)); }

    double llr(int wins, int losses) const {
        const double p1 = win_chance(elo);
        return wins * std::log(2 * p1) + losses * std::log(2 * (1 - p1));
    }
    double upper() const { return std::log((1 - beta) / alpha); }
    double lower() const { return std::log(beta / (1 - alpha)); }

    // +1 when H1 is taken (the first is stronger), -1 when H0 is, 0 while the games must go on.
    int verdict(int wins, int losses) const {
        const double l = llr(wins, losses);
        return l >= upper() ? 1 : l <= lower() ? -1 : 0;
    }
};

} // namespace spatium::statistics
