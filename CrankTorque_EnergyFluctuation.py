import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parameters ----------------
E_number = 395
R1 = (E_number % 6) + 1
R2 = (E_number % 3) + 1
R3 = (E_number % 4) + 1

omega_rpm = 1660       # crank speed in RPM
omega = omega_rpm * 2 * np.pi / 60  # rad/s

Bore = 0.08            # m
Stroke = 0.11          # m
Crank_radius = Stroke / 2
Connection_rod_length = 0.235 # m
Area_of_piston = np.pi * (Bore / 2)**2
Compression_ratio = 8

Mass_of_reciprocating_parts_per_cylinder = 1.5 + R1/10  # kg

# ---------------- Read Data ----------------
df = pd.read_excel("data.xlsx")  # make sure your Excel file has Crank_Angle and R3=4 columns
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
crank_angle_rad = np.radians(crank_angle_deg)

# ---------------- Calculate Mean Torque ----------------
work_per_cycle = np.trapz(TOR, crank_angle_rad)
mean_torque = work_per_cycle / (4 * np.pi)
print(f"Work done per cylinder per cycle: {work_per_cycle:.2f} J")
print(f"Mean torque: {mean_torque:.2f} N·m")

# ---------------- Calculate Energy Fluctuations ----------------
above_mean = TOR - mean_torque
energy_segments = []
labels = []
cum_energy = []
cumulative = 0

sign = np.sign(above_mean[0])
start_idx = 0
segment_idx = 1
for i in range(1, len(TOR)):
    current_sign = np.sign(above_mean[i])
    if current_sign != sign or i == len(TOR)-1:
        segment_energy = np.trapz(above_mean[start_idx:i+1], crank_angle_rad[start_idx:i+1])
        energy_segments.append(segment_energy)
        cumulative += segment_energy
        cum_energy.append(cumulative)
        labels.append(f"a{segment_idx}")
        segment_idx += 1
        start_idx = i
        sign = current_sign

# ---------------- Plot Torque with Shaded Energy ----------------
plt.figure(figsize=(14,7))
plt.plot(crank_angle_deg, TOR, color='blue', lw=2, label="Torque")

# Horizontal and Vertical Lines
plt.axhline(y=0, color='black', lw=1, label="Zero Torque")
plt.axhline(y=mean_torque, color='green', lw=2, linestyle='--', label=f"Mean Torque = {mean_torque:.2f} N·m")

strokes = {"Intake":0, "Compression":180, "Power":360, "Exhaust":540, "Next Intake":720}
for stroke, angle in strokes.items():
    plt.axvline(x=angle, color='red', linestyle='--', lw=1)
    plt.text(angle+5, max(TOR)*0.8, stroke, rotation=90, color='red', verticalalignment='center')

# Shade Energy Areas and Add Labels
start_idx = 0
sign = np.sign(above_mean[0])
segment_idx = 0
y_offset = max(TOR)*0.05  # spacing for labels

for i in range(1, len(TOR)):
    current_sign = np.sign(above_mean[i])
    if current_sign != sign or i == len(TOR)-1:
        x_segment = crank_angle_deg[start_idx:i+1]
        y_segment = TOR[start_idx:i+1]
        mid_angle = x_segment[len(x_segment)//2]

        if energy_segments[segment_idx] >= 0:
            plt.fill_between(x_segment, mean_torque, y_segment, color='orange', alpha=0.3)
            label_y = max(y_segment) + y_offset
            value_y = max(y_segment) + y_offset*1.5
        else:
            plt.fill_between(x_segment, mean_torque, y_segment, color='purple', alpha=0.3)
            label_y = min(y_segment) - y_offset
            value_y = min(y_segment) - y_offset*1.5

        plt.text(mid_angle, label_y, labels[segment_idx], color='black', fontweight='bold',
                 ha='center', va='center', fontsize=10)
        plt.text(mid_angle, value_y, f"{energy_segments[segment_idx]:.2f} J",
                 color='black', ha='center', va='center', fontsize=9)

        start_idx = i
        sign = current_sign
        segment_idx += 1

plt.xlabel("Crank Angle (deg)")
plt.ylabel("Torque (N·m)")
plt.title("Torque vs Crank Angle with Energy Fluctuations")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# ---------------- Stroke-wise Energy and Maximum Fluctuation ----------------
stroke_ranges = {
    "Intake": (0, 180),
    "Compression": (180, 360),
    "Power": (360, 540),
    "Exhaust": (540, 720)
}

stroke_energy = {}
for stroke, (start_deg, end_deg) in stroke_ranges.items():
    idx = np.where((crank_angle_deg >= start_deg) & (crank_angle_deg <= end_deg))[0]
    energy = np.trapz(TOR[idx] - mean_torque, crank_angle_rad[idx])
    stroke_energy[stroke] = energy

stroke_energy_table = pd.DataFrame({
    "Stroke": list(stroke_energy.keys()),
    "Energy (J)": list(stroke_energy.values())
})

max_energy_fluctuation = stroke_energy_table["Energy (J)"].max() - stroke_energy_table["Energy (J)"].min()

print("\nStroke-wise Energy Table:")
print(stroke_energy_table)
print(f"\nMaximum energy fluctuation ΔE: {max_energy_fluctuation:.2f} J")
