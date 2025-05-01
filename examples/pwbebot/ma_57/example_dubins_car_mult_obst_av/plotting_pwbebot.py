import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def read_data(filename):
    df = pd.read_csv(filename)
    return df['Time'].tolist(), df['Value'].tolist()

#── Read trajectory (x₁,x₂) ────────────────────────────────────────────────────────
t_path, x1 = read_data('x1.csv')
_, x2 = read_data('x2.csv')
cp_t_x1, cp_x1 = read_data('x1_controlpoints.csv')
cp_t_x2, cp_x2 = read_data('x2_controlpoints.csv')
cont_t_x1, cont_x1 = read_data('x1_continuity.csv')
cont_t_x2, cont_x2 = read_data('x2_continuity.csv')

#── Read ψ (heading) ─────────────────────────────────────────────────────────────
t_psi, psi = read_data('psi.csv')
cp_t_psi, cp_psi = read_data('psi_controlpoints.csv')
cont_t_psi, cont_psi = read_data('psi_continuity.csv')

#── Read v (speed) ───────────────────────────────────────────────────────────────
t_v, v = read_data('v.csv')
cp_t_v, cp_v = read_data('v_controlpoints.csv')
cont_t_v, cont_v = read_data('v_continuity.csv')

#── Read ω (turn rate) ───────────────────────────────────────────────────────────
t_om, om = read_data('om.csv')
cp_t_om, cp_om = read_data('om_controlpoints.csv')
cont_t_om, cont_om = read_data('om_continuity.csv')

#── Read obstacles ───────────────────────────────────────────────────────────────
obs_df = pd.read_csv('obstacles.csv')
obs_x = obs_df['x'].tolist()
obs_y = obs_df['y'].tolist()
obs_r = obs_df['radius'].tolist()


#── 1) XY path plot ─────────────────────────────────────────────────────────────
plt.figure(figsize=(6,6))
ax = plt.gca()

# trajectory
plt.plot(x1, x2, '-', linewidth=0.8, label='Trajectory')

# control points
plt.scatter(cp_x1, cp_x2, facecolors='none', edgecolors='C0', s=80, label='Control Points')

# continuity points
plt.scatter(cont_x1, cont_x2, facecolors='none', edgecolors='k', s=120, label='Continuity Points')

# obstacles
for xi, yi, ri in zip(obs_x, obs_y, obs_r):
    circle = patches.Circle(
        (xi, yi), ri,
        edgecolor='red', facecolor='none',
        linestyle='--', linewidth=1.2
    )
    ax.add_patch(circle)
# add a single legend entry for obstacles
obstacle_patch = patches.Patch(
    edgecolor='red', facecolor='none',
    linestyle='--', linewidth=1.2, label='Obstacles'
)
ax.legend(handles=ax.get_legend_handles_labels()[0] + [obstacle_patch],
          labels=ax.get_legend_handles_labels()[1] + ['Obstacles'])

plt.xlabel('x₁')
plt.ylabel('x₂')
plt.title('PWBeBOT Path in XY Plane')
plt.axis('equal')
plt.grid(True)


#── 2) ψ (heading) vs time ───────────────────────────────────────────────────────
plt.figure(figsize=(8,4))
plt.plot(t_psi, psi, ':', linewidth=0.8, marker='.', label='ψ')
plt.scatter(cp_t_psi, cp_psi, facecolors='none', edgecolors='C0', s=80, label='Control Points')
plt.scatter(cont_t_psi, cont_psi, facecolors='none', edgecolors='k', s=120, label='Continuity Points')
plt.xlabel('Time')
plt.ylabel('ψ (rad)')
plt.title('Heading ψ over Time')
plt.grid(True)
plt.legend()


#── 3) v (speed) vs time ─────────────────────────────────────────────────────────
plt.figure(figsize=(8,4))
plt.plot(t_v, v, ':', linewidth=0.8, marker='.', label='v')
plt.scatter(cp_t_v, cp_v, facecolors='none', edgecolors='C0', s=80, label='Control Points')
plt.scatter(cont_t_v, cont_v, facecolors='none', edgecolors='k', s=120, label='Continuity Points')
plt.xlabel('Time')
plt.ylabel('v (m/s)')
plt.title('Speed v over Time')
plt.grid(True)
plt.legend()


#── 4) ω (turn rate) vs time ─────────────────────────────────────────────────────
plt.figure(figsize=(8,4))
plt.plot(t_om, om, ':', linewidth=0.8, marker='.', label='ω')
plt.scatter(cp_t_om, cp_om, facecolors='none', edgecolors='C0', s=80, label='Control Points')
plt.scatter(cont_t_om, cont_om, facecolors='none', edgecolors='k', s=120, label='Continuity Points')
plt.xlabel('Time')
plt.ylabel('ω (rad/s)')
plt.title('Angular Rate ω over Time')
plt.grid(True)
plt.legend()


plt.show()

