import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

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

# ---------------- Calculate Mean Torque ----------------
crank_angle_rad = np.radians(crank_angle_deg)
work_per_cycle = np.trapz(TOR, crank_angle_rad)  # Joules
mean_torque = work_per_cycle / (4* np.pi)  # Mean torque (N·m)
print(f"Work done per cylinder per cycle: {work_per_cycle:.2f} J")
print(f"Mean torque: {mean_torque:.2f} N·m")

# ---------------- Plot ----------------
fig, axs = plt.subplots(2,1, figsize=(12,12), sharex=True)

# ----- Top subplot: Full torque curve with work area -----
axs[0].fill_between(crank_angle_deg, TOR, 0, color='skyblue', alpha=0.3, label="Work Area")
axs[0].plot(crank_angle_deg, TOR, color='blue', lw=2, label="Torque")
axs[0].axhline(y=mean_torque, color='green', linestyle='--', lw=2, label=f"Mean Torque = {mean_torque:.2f} N·m")
axs[0].axhline(y=0, color='black', lw=1)

# Vertical stroke lines
strokes = {"Intake":0, "Compression":180, "Power":360, "Exhaust":540, "Next Intake":720}
for stroke, angle in strokes.items():
    axs[0].axvline(x=angle, color='red', linestyle='--', lw=1)
    axs[0].text(angle+5, max(TOR)*0.8, stroke, rotation=90, color='red', verticalalignment='center')

axs[0].set_ylabel("Torque (N·m)")
axs[0].set_title("Torque vs Crank Angle with Mean Torque")
axs[0].grid(True)
axs[0].legend()

# ----- Bottom subplot: Areas above and below mean torque with labels -----
axs[1].fill_between(crank_angle_deg, TOR, mean_torque, where=(TOR>mean_torque), color='lightgreen', alpha=0.5, label="Above Mean Torque")
axs[1].fill_between(crank_angle_deg, TOR, mean_torque, where=(TOR<mean_torque), color='salmon', alpha=0.5, label="Below Mean Torque")
axs[1].plot(crank_angle_deg, TOR, color='blue', lw=2, label="Torque")
axs[1].axhline(y=mean_torque, color='green', linestyle='--', lw=2, label=f"Mean Torque = {mean_torque:.2f} N·m")
axs[1].axhline(y=0, color='black', lw=1)

# Vertical stroke lines
for stroke, angle in strokes.items():
    axs[1].axvline(x=angle, color='red', linestyle='--', lw=1)

# ---- Label continuous segments A1, A2,... and B1, B2,... with their area ----
above = TOR > mean_torque
below = TOR < mean_torque

def label_segments_with_area(condition, torque_array, angle_array, label_prefix):
    labels = []
    start_idx = None
    seg_num = 1
    areas = []
    for i in range(len(condition)):
        if condition[i] and start_idx is None:
            start_idx = i
        elif not condition[i] and start_idx is not None:
            end_idx = i
            # Midpoint for label
            mid_idx = (start_idx + end_idx)//2
            # Numerical area using trapezoidal rule
            area = np.trapz(torque_array[start_idx:end_idx] - mean_torque, np.radians(angle_array[start_idx:end_idx]))
            areas.append(area)
            labels.append((mid_idx, f"{label_prefix}{seg_num}", area))
            seg_num += 1
            start_idx = None
    if start_idx is not None:
        end_idx = len(condition)-1
        mid_idx = (start_idx + end_idx)//2
        area = np.trapz(torque_array[start_idx:end_idx] - mean_torque, np.radians(angle_array[start_idx:end_idx]))
        areas.append(area)
        labels.append((mid_idx, f"{label_prefix}{seg_num}", area))
    return labels

# Label above mean
for idx, lab, area in label_segments_with_area(above, TOR, crank_angle_deg, "A"):
    axs[1].text(crank_angle_deg[idx], TOR[idx]+max(TOR)*0.02, f"{lab}\n{area:.1f}", color='darkgreen', ha='center')

# Label below mean
for idx, lab, area in label_segments_with_area(below, TOR, crank_angle_deg, "B"):
    axs[1].text(crank_angle_deg[idx], TOR[idx]-max(TOR)*0.05, f"{lab}\n{area:.1f}", color='darkred', ha='center')

axs[1].set_xlabel("Crank Angle (deg)")
axs[1].set_ylabel("Torque (N·m)")
axs[1].set_title("Torque Areas Above and Below Mean Torque with Numerical Values")
axs[1].grid(True)
axs[1].legend()
axs[1].set_xlim(0,720)

plt.tight_layout()
plt.show()
