import pandas as pd
import matplotlib.pyplot as plt

def read_data(filename):
    df = pd.read_csv(filename)
    return df['Time'].tolist(), df['Value'].tolist()

# Read data for z, w, theta, q
times_z, values_z = read_data('z.csv')
times_w, values_w = read_data('w.csv')
times_theta, values_theta = read_data('theta.csv')
times_q, values_q = read_data('q.csv')

cp_times_z, cp_values_z = read_data('z_controlpoints.csv')
cp_times_w, cp_values_w = read_data('w_controlpoints.csv')
cp_times_theta, cp_values_theta = read_data('theta_controlpoints.csv')
cp_times_q, cp_values_q = read_data('q_controlpoints.csv')

# Read data for delta_v, delta_s, delta_m
times_delta_v, values_delta_v = read_data('delta_v.csv')
times_delta_s, values_delta_s = read_data('delta_s.csv')
times_delta_m, values_delta_m = read_data('delta_m.csv')

cp_times_delta_v, cp_values_delta_v = read_data('delta_v_controlpoints.csv')
cp_times_delta_s, cp_values_delta_s = read_data('delta_s_controlpoints.csv')
cp_times_delta_m, cp_values_delta_m = read_data('delta_m_controlpoints.csv')

# Plotting z, w, theta, q
plt.figure(figsize=(10, 6))

line_width = 0.5

plt.plot(times_z, values_z, marker='.', linestyle=':', color='blue', linewidth=line_width, label='z')
plt.plot(times_w, values_w, marker='.', linestyle=':', color='red', linewidth=line_width, label='w')
plt.plot(times_theta, values_theta, marker='.', linestyle=':', color='green', linewidth=line_width, label='theta')
plt.plot(times_q, values_q, marker='.', linestyle=':', color='purple', linewidth=line_width, label='q')

plt.scatter(cp_times_z, cp_values_z, color='white', s=70, label='z Control Points', edgecolor='blue') 
plt.scatter(cp_times_w, cp_values_w, color='white', s=70, label='w Control Points', edgecolor='red') 
plt.scatter(cp_times_theta, cp_values_theta, color='white', s=70, label='theta Control Points', edgecolor='green') 
plt.scatter(cp_times_q, cp_values_q, color='white', s=70, label='q Control Points', edgecolor='purple') 

plt.title('BeBOT Results: z, w, theta, q')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()

# Plotting delta_v, delta_s, delta_m
plt.figure(figsize=(10, 6))

plt.plot(times_delta_v, values_delta_v, marker='.', linestyle=':', color='blue', linewidth=line_width, label='delta_v')
plt.plot(times_delta_s, values_delta_s, marker='.', linestyle=':', color='red', linewidth=line_width, label='delta_s')
plt.plot(times_delta_m, values_delta_m, marker='.', linestyle=':', color='green', linewidth=line_width, label='delta_m')

plt.scatter(cp_times_delta_v, cp_values_delta_v, color='white', s=70, label='delta_v Control Points', edgecolor='blue') 
plt.scatter(cp_times_delta_s, cp_values_delta_s, color='white', s=70, label='delta_s Control Points', edgecolor='red') 
plt.scatter(cp_times_delta_m, cp_values_delta_m, color='white', s=70, label='delta_m Control Points', edgecolor='green') 

plt.title('BeBOT Results: delta_v, delta_s, delta_m')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()

# Show both plots at the end
plt.show()

