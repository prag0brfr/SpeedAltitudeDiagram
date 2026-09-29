#!/usr/bin/env python3
"""Interactive ISA altitude–TAS diagram with iso-CAS and iso-Mach curves.

Run: python altitude_tas_diagram.py
Click inside the plot to read altitude, TAS, CAS, and Mach at that point.
Requires: numpy and matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt

# Plot defaults (edit these to suit the flight envelope of interest).
MAX_ALTITUDE_FT = 45_000
MAX_TAS_KT = 500
CAS_CONTOURS_KT = np.arange(100, 501, 25)
MACH_CONTOURS = np.arange(0.2, 1.001, 0.05)

# ISA constants
G0 = 9.80665                 # m/s²
R = 287.05287                # J/(kg·K)
GAMMA = 1.4
T0 = 288.15                  # K
P0 = 101325.0                # Pa
LAPSE = 0.0065               # K/m, tropospheric lapse rate
TROPOPAUSE_M = 11000.0       # m
T11 = T0 - LAPSE * TROPOPAUSE_M
P11 = P0 * (T11 / T0) ** (G0 / (R * LAPSE))
A0 = np.sqrt(GAMMA * R * T0)
MPS_TO_KT = 1.9438444924406
FT_TO_M = 0.3048


def isa_atmosphere(altitude_m):
    """Return ISA static temperature (K) and pressure (Pa), 0–20 km."""
    h = np.asarray(altitude_m, dtype=float)
    if np.any((h < 0) | (h > 20000)):
        raise ValueError("ISA altitude must be between 0 and 20,000 m.")
    t = np.where(h <= TROPOPAUSE_M, T0 - LAPSE * h, T11)
    p_trop = P0 * (t / T0) ** (G0 / (R * LAPSE))
    p_strat = P11 * np.exp(-G0 * (h - TROPOPAUSE_M) / (R * T11))
    p = np.where(h <= TROPOPAUSE_M, p_trop, p_strat)
    return t, p


def impact_pressure(mach, static_pressure):
    """Compressible-flow impact pressure qc for subsonic Mach."""
    return static_pressure * ((1.0 + (GAMMA - 1.0) * mach**2 / 2.0) ** (GAMMA / (GAMMA - 1.0)) - 1.0)


def cas_from_mach_pressure(mach, static_pressure):
    """Convert Mach and static pressure to calibrated airspeed in knots."""
    qc = impact_pressure(mach, static_pressure)
    m_cas_sq = (2.0 / (GAMMA - 1.0)) * (
        (qc / P0 + 1.0) ** ((GAMMA - 1.0) / GAMMA) - 1.0
    )
    return A0 * np.sqrt(np.maximum(m_cas_sq, 0.0)) * MPS_TO_KT


def mach_from_tas( tas_kt, temperature_k):
    """Convert true airspeed in knots to Mach."""
    return (np.asarray(tas_kt) / MPS_TO_KT) / np.sqrt(GAMMA * R * temperature_k)


def main():
    max_alt_m = MAX_ALTITUDE_FT * FT_TO_M
    if max_alt_m > 20000:
        raise ValueError("MAX_ALTITUDE_FT must not exceed about 65,617 ft (20,000 m).")

    # Grid for the altitude-dependent performance curves.
    alt_m = np.linspace(0.0, max_alt_m, 500)
    temp_k, pressure_pa = isa_atmosphere(alt_m)
    sound_speed_kt = np.sqrt(GAMMA * R * temp_k) * MPS_TO_KT
    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)

    # Iso-Mach lines: TAS = Mach × local speed of sound.
    for mach in MACH_CONTOURS:
        tas_kt = mach * sound_speed_kt
        visible = tas_kt <= MAX_TAS_KT
        if np.any(visible):
            major = np.isclose(mach * 10, round(mach * 10))
            ax.plot(tas_kt[visible], alt_m[visible] / FT_TO_M,
                    color="tab:orange", lw=1.2 if major else 0.7,
                    alpha=0.9 if major else 0.5, ls="--")
            if major:
                visible_idx = np.flatnonzero(visible)
                # Place Mach labels inside the plot, away from the title.
                idx = visible_idx[len(visible_idx) // 2]
                ax.annotate(f"M {mach:.1f}", (tas_kt[idx], alt_m[idx] / FT_TO_M),
                            xytext=(4, 3), textcoords="offset points", fontsize=8,
                            color="darkorange")

    # Iso-CAS lines: at each altitude, solve the compressible impact-pressure relation.
    for cas_kt in CAS_CONTOURS_KT:
        cas_mps = cas_kt / MPS_TO_KT
        qc = P0 * ((1.0 + (GAMMA - 1.0) * (cas_mps / A0) ** 2 / 2.0)
                   ** (GAMMA / (GAMMA - 1.0)) - 1.0)
        mach_sq = (2.0 / (GAMMA - 1.0)) * (
            (qc / pressure_pa + 1.0) ** ((GAMMA - 1.0) / GAMMA) - 1.0
        )
        tas_kt = np.sqrt(np.maximum(mach_sq, 0.0)) * sound_speed_kt
        visible = tas_kt <= MAX_TAS_KT
        if np.any(visible):
            major = cas_kt % 50 == 0
            ax.plot(tas_kt[visible], alt_m[visible] / FT_TO_M,
                    color="tab:blue", lw=1.15 if major else 0.7,
                    alpha=0.9 if major else 0.5)
            if major:
                # Keep major CAS labels aligned at 20,000 ft (FL200).
                fl200_ft = min(20_000.0, MAX_ALTITUDE_FT)
                fl200_m = fl200_ft * FT_TO_M
                tas_at_fl200 = float(np.interp(fl200_m, alt_m, tas_kt))
                if tas_at_fl200 <= MAX_TAS_KT:
                    ax.annotate(f"{cas_kt:.0f} CAS",
                                (tas_at_fl200, fl200_ft),
                                xytext=(4, -9), textcoords="offset points", fontsize=8,
                                color="tab:blue")

    ax.set_xlim(0, MAX_TAS_KT)
    ax.set_ylim(0, MAX_ALTITUDE_FT)
    ax.set_xlabel("True airspeed (kt TAS)")
    ax.set_ylabel("Pressure altitude (ft, ISA)")
    ax.set_title("Altitude–TAS Diagram — ISA, click to read the flight condition")
    ax.grid(True, alpha=0.28)

    # Persistent selected-point marker and readout.
    marker, = ax.plot([], [], "ro", ms=7, zorder=10)
    readout = ax.text(
        0.99, 0.02, "Click inside the plot to select a point.",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.92,
                  edgecolor="0.6"),
    )

    def on_click(event):
        if event.inaxes is not ax or event.xdata is None or event.ydata is None:
            return
        tas_kt = float(np.clip(event.xdata, 0.0, MAX_TAS_KT))
        altitude_ft = float(np.clip(event.ydata, 0.0, MAX_ALTITUDE_FT))
        altitude_m = altitude_ft * FT_TO_M
        temp, pressure = isa_atmosphere(altitude_m)
        mach = float(mach_from_tas(tas_kt, temp))
        cas_kt = float(cas_from_mach_pressure(mach, pressure))
        marker.set_data([tas_kt], [altitude_ft])
        readout.set_text(
            f"Altitude: {altitude_ft:,.0f} ft\n"
            f"TAS: {tas_kt:.1f} kt\n"
            f"CAS: {cas_kt:.1f} kt\n"
            f"Mach: {mach:.3f}"
        )
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("button_press_event", on_click)
    plt.show()


if __name__ == "__main__":
    main()
