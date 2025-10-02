# Blower model (1-D, steady, isothermal, compressible air)
# Now also computes electrical input power from the fan curve

import math
import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# User-settable parameters
# -----------------------
L     = 0.5          # blower/duct length [m]
D     = 0.10         # hydraulic diameter [m]
A     = math.pi*(D**2)/4.0  # cross-sectional area [m^2]
f_D   = 0.02         # Darcy friction factor [-] (assumed constant)
T     = 298.15       # temperature [K] (isothermal)
R     = 287.0        # specific gas constant for dry air [J/(kg·K)]
p_in  = 101325.0     # inlet static pressure [Pa] (near 1 atm)

# ---- Fan (blower) characteristic: total pressure rise vs inlet volumetric flow
dp0       = 800.0     # shut-off pressure rise [Pa] at Q=0
Q_max     = 0.25      # zero-pressure-flow (max flow) [m^3/s]
def dp_fan_of_Qin(Qin):
    """Fan total pressure rise (Pa) as a function of inlet volumetric flow (m^3/s)."""
    x = 1.0 - (Qin/Q_max)**2
    return dp0*max(0.0, x)

# ---- Hardware/drive efficiencies (fractions 0-1)
eta_fan        = 0.65   # aerodynamic (total-to-total) efficiency of the fan
eta_drive      = 0.97   # belt/gear/direct drive (1.0 for direct-coupled)
eta_motor      = 0.93   # motor efficiency
eta_controller = 0.98   # VFD/speed controller (1.0 if none)

# Integration/grid
N  = 400
z  = np.linspace(0.0, L, N+1)
dz = z[1]-z[0]

# -----------------------------------------
# Core integrator: given m_dot, integrate p(z)
# -----------------------------------------
def integrate_pressure_profile(m_dot):
    """
    Integrates dp/dz = -(f_D/D) * rho * u^2 + S_fan
      with rho = p/(R T), u = m_dot/(rho A), S_fan = Δp_fan(Q_in)/L
    Returns arrays p(z), Q(z), and the fan total rise used.
    """
    rho_in = p_in/(R*T)
    u_in   = m_dot/(rho_in*A)
    Q_in   = A*u_in

    dp_fan_total = dp_fan_of_Qin(Q_in)  # total rise Pa from fan characteristic
    S_fan = dp_fan_total / L            # distributed source [Pa/m]

    p = np.empty_like(z)
    p[0] = p_in
    for i in range(N):
        rho = p[i]/(R*T)
        u   = m_dot/(rho*A)
        dp_dz = - (f_D/D) * rho * u*u + S_fan
        p[i+1] = p[i] + dp_dz*dz

    rho_z = p/(R*T)
    Q_z   = m_dot / rho_z
    return p, Q_z, dp_fan_total, Q_in

# ----------------------------------------------------------
# 1) Build a p–Q curve and compute power over a sweep of Qin
# ----------------------------------------------------------
Qin_sweep = np.linspace(0.01*Q_max, 0.95*Q_max, 60)
dP_net    = []     # p_out - p_in
P_shaft   = []     # air/shaft power = Qin * dp_fan_total  [W]
P_elec    = []     # electrical input power                 [W]

for Qin in Qin_sweep:
    rho_in = p_in/(R*T)
    m_dot  = rho_in * Qin
    p_prof, Q_prof, dp_fan_total, Q_in_used = integrate_pressure_profile(m_dot)

    # Net pressure change across the length (includes friction): p_out - p_in
    dP_net.append(p_prof[-1] - p_prof[0])

    # Power from fan curve at inlet condition (wire-to-air uses Qin at inlet)
    P_air   = Q_in_used * dp_fan_total                # useful power to the air [W]
    Pshaft  = P_air / max(eta_fan, 1e-9)             # shaft input to fan [W]
    Pelec   = Pshaft / (max(eta_drive*eta_motor*eta_controller, 1e-9))  # [W]

    P_shaft.append(Pshaft)
    P_elec.append(Pelec)

dP_net  = np.array(dP_net)
P_shaft = np.array(P_shaft)
P_elec  = np.array(P_elec)

# -------------------------------------------------------
# 2) Choose an operating point and plot Q vs z and p vs z
# -------------------------------------------------------
# Choose mid-range operating point (you can change this)
Qin_op   = 0.5*(Qin_sweep[0] + Qin_sweep[-1])
rho_in   = p_in/(R*T)
m_dot_op = rho_in * Qin_op
p_op, Qz_op, dp_fan_op, Qin_used_op = integrate_pressure_profile(m_dot_op)

# Compute powers at operating point
P_air_op   = Qin_used_op * dp_fan_op
P_shaft_op = P_air_op / max(eta_fan, 1e-9)
P_elec_op  = P_shaft_op / (max(eta_drive*eta_motor*eta_controller, 1e-9))

# -----------------
# Plotting section
# -----------------
# (1) Blower p–Q curve (net across length)
plt.figure()
plt.plot(Qin_sweep, dP_net/1000.0, lw=2)
plt.axhline(0.0, ls='--')
plt.xlabel('Inlet volumetric flow $Q_{in}$ [m$^3$/s]')
plt.ylabel('Net pressure change $\\Delta p = p_{out}-p_{in}$ [kPa]')
plt.title('Blower p–Q curve (net across length)')
plt.grid(True)

# (2) Electrical input power vs flow (computed from fan curve and efficiencies)
plt.figure()
plt.plot(Qin_sweep, P_elec/1000.0, lw=2)
plt.xlabel('Inlet volumetric flow $Q_{in}$ [m$^3$/s]')
plt.ylabel('Electrical input power $P_{elec}$ [kW]')
plt.title('Electrical input power vs inlet flow')
plt.grid(True)

# (3) Pressure profile at operating point
plt.figure()
plt.plot(z, p_op/1000.0, lw=2)
plt.xlabel('Axial position z [m]')
plt.ylabel('Static pressure p(z) [kPa]')
plt.title('Pressure profile along blower (operating point)')
plt.grid(True)

# (4) Volumetric flow profile at operating point (compressibility)
plt.figure()
plt.plot(z, Qz_op, lw=2)
plt.xlabel('Axial position z [m]')
plt.ylabel('Volumetric flow $Q(z)$ [m$^3$/s]')
plt.title('Volumetric flow along blower (operating point)')
plt.grid(True)

# (5) Bar chart: shaft vs electrical power at operating point
plt.figure()
labels = ['Air power', 'Shaft power', 'Electrical power']
vals   = [P_air_op/1000.0, P_shaft_op/1000.0, P_elec_op/1000.0]
plt.bar(labels, vals)
plt.ylabel('Power [kW]')
plt.title('Power breakdown at operating point')
plt.grid(axis='y')

plt.show()

# -----------------
# Print a quick numeric summary for the chosen operating point
# -----------------
print("\n--- Operating point summary ---")
print(f"Q_in used          : {Qin_used_op:.4f} m^3/s")
print(f"Fan Δp_total       : {dp_fan_op:.1f} Pa")
print(f"Air power (Q*Δp)   : {P_air_op/1000.0:.3f} kW")
print(f"Shaft efficiency   : eta_fan = {eta_fan:.3f}")
print(f"Drive*Motor*Ctrl   : {eta_drive:.3f} * {eta_motor:.3f} * {eta_controller:.3f} = {(eta_drive*eta_motor*eta_controller):.3f}")
print(f"Shaft power        : {P_shaft_op/1000.0:.3f} kW")
print(f"Electrical power   : {P_elec_op/1000.0:.3f} kW")
