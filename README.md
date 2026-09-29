# Interactive Altitude–TAS Diagram

A small interactive Python tool for exploring the relationship between pressure altitude, true airspeed (TAS), calibrated airspeed (CAS), and Mach number under International Standard Atmosphere (ISA) conditions.

The plot shows iso-CAS and iso-Mach curves. Click anywhere inside the plot to display the altitude, TAS, CAS, and Mach number for that point.

## Requirements

- Python 3.9 or later
- NumPy
- Matplotlib

## Install

Install the dependencies with pip:

```bash
python -m pip install numpy matplotlib
```

## Run

Place `README.md` alongside `SpeedAltitudeDiagram.py`, then run:

```bash
python SpeedAltitudeDiagram.py
```

Click within the plot to select a flight condition. The selected point is marked in red, and its four values appear in the lower-right corner. Click another point to update the readout. Close the plot window to exit.

## Plot settings

Edit the constants near the top of `SpeedAltitudeDiagram.py` to adjust the plotted range and reference curves:

```python
MAX_ALTITUDE_FT = 45_000
MAX_TAS_KT = 500
CAS_CONTOURS_KT = np.arange(100, 501, 25)
MACH_CONTOURS = np.arange(0.2, 1.001, 0.05)
```

The current defaults cover 0–45,000 ft and 0–500 kt TAS, with CAS contours every 25 kt and Mach contours every 0.05. Major reference curves are emphasized and labeled; there is no legend. CAS labels are positioned at FL200 (or at the maximum plotted altitude if it is below 20,000 ft), and Mach labels are placed within the plot.

`MAX_ALTITUDE_FT` must not exceed approximately 65,617 ft (20,000 m), the upper limit of the atmosphere model used by the script. Curves outside the selected TAS range are clipped from view.

## Assumptions and limitations

- The atmosphere is modeled using ISA from sea level through 20,000 m: a linear temperature lapse rate up to 11,000 m, followed by an isothermal layer.
- CAS is calculated using the compressible-flow impact-pressure relation, referenced to standard sea-level pressure and temperature.
- The equations and curve generation assume subsonic flow. Treat results near or above Mach 1 as indicative, not as validated supersonic air-data calculations.
- The displayed values are idealized atmosphere/airspeed conversions, not aircraft performance limits or flight guidance. Instrument errors, position error, non-ISA weather, and aircraft-specific effects are not modeled.

## License

No license is specified. Add a `LICENSE` file before redistributing this repository if you want to grant others explicit reuse rights.
