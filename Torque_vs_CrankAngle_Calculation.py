
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parameters ----------------
E_number = 395
R1 = (E_number % 6) + 1
R2 = (E_number % 3) + 1
R3 = (E_number % 4) + 1

omega_rpm = 2000       # crank speed in RPM
omega = omega_rpm * 2 * np.pi / 60  # rad/s

Bore = 0.08            # m
Stroke = 0.11          # m
Crank_radius = Stroke / 2
Connection_rod_length = 0.235 # m
Area_of_piston = np.pi * (Bore / 2)**2
Compression_ratio = 8

Mass_of_reciprocating_parts_per_cylinder = 1.5 + R1/10  # kg

# ---------------- Read Data ----------------
df = pd.read_excel("data.xlsx")
crank_angle_deg = np.array(df["Crank_Angle"])
pressure_bar = np.array(df["R3 = 4"])
pressure_Pa = pressure_bar * 1e5  # convert bar to Pa

# ---------------- Functions ----------------
def Q(theta_rad):
    return Mass_of_reciprocating_parts_per_cylinder * omega**2 * Crank_radius * (
        np.cos(theta_rad) + (np.cos(2*theta_rad)/Compression_ratio)
    )

def Torque(p, theta_rad):
    sin_theta = np.sin(theta_rad)
    term = 2 * Crank_radius * np.sqrt(Connection_rod_length**2 - (Crank_radius * sin_theta)**2)
    torque = (p * Area_of_piston - Q(theta_rad)) * Crank_radius * (
        sin_theta + (Crank_radius * sin_theta) / term
    )
    return torque

# ---------------- Calculate Torque ----------------
TOR = np.array([Torque(p, np.radians(theta)) for p, theta in zip(pressure_Pa, crank_angle_deg)])

# ---------------- Plot ----------------
plt.figure(figsize=(12,6))
plt.plot(crank_angle_deg, TOR, lw=2, label="Torque")

# Zero torque line
plt.axhline(y=0, lw=1, label="Zero Torque")

# Vertical lines for strokes
strokes = {"Intake":0, "Compression":180, "Power":360, "Exhaust":540, "Next Intake":720}
for stroke, angle in strokes.items():
    plt.axvline(x=angle, linestyle='--', lw=1)
    plt.text(angle+5, max(TOR)*0.8, stroke,
             rotation=90, verticalalignment='center')

# ---------------- Vertical line at 570° ----------------
angle_570 = 570

# Interpolate torque value at 570°
torque_570 = np.interp(angle_570, crank_angle_deg, TOR)

# Plot vertical line
plt.axvline(
    x=angle_570,
    linestyle='-.',
    lw=2,
    label=f"θ = 570°, T = {torque_570:.2f} N·m"
)

# Optional marker at intersection
plt.plot(angle_570, torque_570, 'o')

# Labels and title
plt.xlabel("Crank Angle (deg)")
plt.ylabel("Torque (N·m)")
plt.title("Torque vs Crank Angle")
plt.grid(True)

# Legend
plt.legend()

plt.show()
