import pandas as pd
import matplotlib.pyplot as plt

def read_data(filename):
    df = pd.read_csv(filename)
    return df['Time'].tolist(), df['Value'].tolist()

# 1) Read data for states z, theta, w, q and their control points
times_z, values_z         = read_data('z.csv')
cp_times_z, cp_values_z   = read_data('z_controlpoints.csv')

times_theta, values_theta = read_data('theta.csv')
cp_times_theta, cp_values_theta = read_data('theta_controlpoints.csv')

times_w, values_w         = read_data('w.csv')
cp_times_w, cp_values_w   = read_data('w_controlpoints.csv')

times_q, values_q         = read_data('q.csv')
cp_times_q, cp_values_q   = read_data('q_controlpoints.csv')

# 2) Read data for states y, psi, v, r and their control points
times_y, values_y         = read_data('y.csv')
cp_times_y, cp_values_y   = read_data('y_controlpoints.csv')

times_psi, values_psi     = read_data('psi.csv')
cp_times_psi, cp_values_psi = read_data('psi_controlpoints.csv')

times_v, values_v         = read_data('u.csv')
cp_times_v, cp_values_v   = read_data('u_controlpoints.csv')

times_r, values_r         = read_data('r.csv')
cp_times_r, cp_values_r   = read_data('r_controlpoints.csv')

# 3) Read data for states x, y (y already read above) and their control points
times_x, values_x         = read_data('x.csv')
cp_times_x, cp_values_x   = read_data('x_controlpoints.csv')
# (y and its control points: times_y, values_y, cp_times_y, cp_values_y)

# 4) Read data for control inputs delta_v, delta_s, delta_m, delta_h, delta_n and their control points
times_delta_v, values_delta_v = read_data('delta_v.csv')
cp_times_delta_v, cp_values_delta_v = read_data('delta_v_controlpoints.csv')

times_delta_s, values_delta_s = read_data('delta_s.csv')
cp_times_delta_s, cp_values_delta_s = read_data('delta_s_controlpoints.csv')

times_delta_m, values_delta_m = read_data('delta_m.csv')
cp_times_delta_m, cp_values_delta_m = read_data('delta_m_controlpoints.csv')

times_delta_h, values_delta_h = read_data('delta_h.csv')
cp_times_delta_h, cp_values_delta_h = read_data('delta_h_controlpoints.csv')

times_delta_n, values_delta_n = read_data('delta_n.csv')
cp_times_delta_n, cp_values_delta_n = read_data('delta_n_controlpoints.csv')

line_width = 0.5

# —— Plot 1: z, theta, w, q and their control points —— 
plt.figure(figsize=(10, 6))
plt.plot(times_z, values_z,         marker='.', linestyle=':', linewidth=line_width, label='z')
plt.plot(times_theta, values_theta, marker='.', linestyle=':', linewidth=line_width, label='theta')
plt.plot(times_w, values_w,         marker='.', linestyle=':', linewidth=line_width, label='w')
plt.plot(times_q, values_q,         marker='.', linestyle=':', linewidth=line_width, label='q')

plt.scatter(cp_times_z,   cp_values_z,   s=70, label='z Control Points')
plt.scatter(cp_times_theta, cp_values_theta, s=70, label='theta Control Points')
plt.scatter(cp_times_w,   cp_values_w,   s=70, label='w Control Points')
plt.scatter(cp_times_q,   cp_values_q,   s=70, label='q Control Points')

plt.title('BeBOT States: z, theta, w, q')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()
plt.tight_layout()

# —— Plot 2: y, psi, v, r and their control points —— 
plt.figure(figsize=(10, 6))
plt.plot(times_psi, values_psi,     marker='.', linestyle=':', linewidth=line_width, label='psi')
plt.plot(times_v, values_v,         marker='.', linestyle=':', linewidth=line_width, label='u')
plt.plot(times_r, values_r,         marker='.', linestyle=':', linewidth=line_width, label='r')


plt.scatter(cp_times_v,   cp_values_v,   s=70, label='u Control Points')
plt.scatter(cp_times_psi, cp_values_psi, s=70, label='psi Control Points')
plt.scatter(cp_times_r,   cp_values_r,   s=70, label='r Control Points')

plt.title('BeBOT States: y, psi, v, r')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()
plt.tight_layout()

# —— Plot 3: x, y and their control points —— 
plt.figure(figsize=(10, 6))
plt.plot(times_x, values_x, marker='.', linestyle=':', linewidth=line_width, label='x')
plt.plot(times_y, values_y, marker='.', linestyle=':', linewidth=line_width, label='y')

plt.scatter(cp_times_x, cp_values_x, s=70, label='x Control Points')
plt.scatter(cp_times_y, cp_values_y, s=70, label='y Control Points')

plt.title('BeBOT States: x, y')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()
plt.tight_layout()

# —— Plot 4: delta_v, delta_s, delta_m, delta_h, delta_n and their control points —— 
plt.figure(figsize=(10, 6))
plt.plot(times_delta_v, values_delta_v, marker='.', linestyle=':', linewidth=line_width, label='delta_v')
plt.plot(times_delta_s, values_delta_s, marker='.', linestyle=':', linewidth=line_width, label='delta_s')
plt.plot(times_delta_m, values_delta_m, marker='.', linestyle=':', linewidth=line_width, label='delta_m')
plt.plot(times_delta_h, values_delta_h, marker='.', linestyle=':', linewidth=line_width, label='delta_h')
plt.plot(times_delta_n, values_delta_n, marker='.', linestyle=':', linewidth=line_width, label='delta_n')

plt.scatter(cp_times_delta_v, cp_values_delta_v, s=70, label='δv Control Points')
plt.scatter(cp_times_delta_s, cp_values_delta_s, s=70, label='δs Control Points')
plt.scatter(cp_times_delta_m, cp_values_delta_m, s=70, label='δm Control Points')
plt.scatter(cp_times_delta_h, cp_values_delta_h, s=70, label='δh Control Points')
plt.scatter(cp_times_delta_n, cp_values_delta_n, s=70, label='δn Control Points')

plt.title('BeBOT Control Inputs: δv, δs, δm, δh, δn')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()

