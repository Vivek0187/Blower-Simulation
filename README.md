# Blower-Simulation

# Blower Model (1-D, Isothermal, Compressible Air)

This project simulates a simple **blower/duct system** using a 1-D steady isothermal compressible model.  
It combines **continuity** and **momentum** equations with **Darcy–Weisbach friction** and a distributed fan source to evaluate blower performance.

---

## Features
- **Pressure–Flow (p–Q) Curve:** Net pressure change across the blower vs. inlet volumetric flow.
- **Flow Distribution:** Volumetric flow rate variation along the duct/blower length.
- **Pressure Profile:** Static pressure variation along the axial coordinate.
- **Power Calculations:**
  - Air (useful) power
  - Shaft power (accounts for fan efficiency)
  - Electrical input power (wire-to-air, accounts for drive, motor, and controller efficiency)
- **Plots:**
  - p–Q curve  
  - Electrical input power vs. flow  
  - Pressure profile along blower  
  - Flow profile along blower  
  - Power breakdown (Air vs. Shaft vs. Electrical)

---

## File Description
- `blower_model.py` → Main Python script implementing the model.
- `README.md` → This file.

---

## User Inputs
You can adjust the following parameters in the script:

- **Geometry & Flow**
  - `L` → duct length [m]  
  - `D` → hydraulic diameter [m]  
  - `f_D` → Darcy friction factor [-]

- **Air Properties**
  - `T` → isothermal temperature [K]  
  - `R` → gas constant for dry air [J/kg·K]  
  - `p_in` → inlet static pressure [Pa]  

- **Fan Characteristic**
  - `dp0` → shut-off pressure rise [Pa]  
  - `Q_max` → max flow at zero pressure [m³/s]

- **Efficiencies**
  - `eta_fan` → aerodynamic fan efficiency  
  - `eta_drive` → transmission (belt/gear) efficiency  
  - `eta_motor` → motor efficiency  
  - `eta_controller` → speed controller efficiency  

---

## Outputs
1. **p–Q curve**: Shows operating points where net Δp = 0.  
2. **Electrical input power vs flow**: Energy requirement for different flows.  
3. **Pressure profile**: Static pressure variation along duct length.  
4. **Flow profile**: Compressibility effects on volumetric flow along the blower.  
5. **Bar chart**: Air, shaft, and electrical power comparison at chosen operating point.  

At runtime, a **numeric summary** is also printed for the selected operating point (Q, Δp, Air Power, Shaft Power, Electrical Power).

---

--- Operating point summary ---
Q_in used          : 0.1184 m^3/s
Fan Δp_total       : 600.5 Pa
Air power (Q*Δp)   : 0.071 kW
Shaft efficiency   : eta_fan = 0.650
Drive*Motor*Ctrl   : 0.970 * 0.930 * 0.980 = 0.885
Shaft power        : 0.110 kW
Electrical power   : 0.124 kW


## 📊 Example Plots

### Pressure–Flow (p–Q) Curve
![p–Q curve](pq_curve.png)

### Electrical Input Power vs Flow
![Power curve](power_vs_flow.png)

### Pressure Profile Along Blower
![Pressure profile](pressure_profile.png)

### Flow Distribution Along Blower
![Flow profile](flow_profile.png)

### Power Breakdown
![Power breakdown](power_breakdown.png)


##  How to Run
```bash
pip install numpy matplotlib
python blower_model.py
```
## Author-Vivek Shabareesh

