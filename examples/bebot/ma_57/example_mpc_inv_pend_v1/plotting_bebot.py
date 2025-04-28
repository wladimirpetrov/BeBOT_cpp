import pandas as pd
import matplotlib.pyplot as plt

def read_data(filename):
    df = pd.read_csv(filename)
    return df['Time'].tolist(), df['Value'].tolist()

# Read data for z, w, theta, q
times_theta, values_theta = read_data('theta.csv')
times_thetadot, values_thetadot = read_data('thetadot.csv')

cp_times_theta, cp_values_theta = read_data('theta_controlpoints.csv')
cp_times_thetadot, cp_values_thetadot = read_data('thetadot_controlpoints.csv')

# Read data for delta_v, delta_s, delta_m
times_f, values_f = read_data('f.csv')

cp_times_f, cp_values_f = read_data('f_controlpoints.csv')

# Plotting z, w, theta, q
plt.figure(figsize=(10, 6))

line_width = 0.5

plt.plot(times_theta, values_theta, marker='.', linestyle=':', color='blue', linewidth=line_width, label='theta')
plt.plot(times_thetadot, values_thetadot, marker='.', linestyle=':', color='red', linewidth=line_width, label='thetadot')

plt.scatter(cp_times_theta, cp_values_theta, color='white', s=70, label='theta Control Points', edgecolor='blue') 
plt.scatter(cp_times_thetadot, cp_values_thetadot, color='white', s=70, label='thetadot Control Points', edgecolor='red') 

plt.title('BeBOT Results: theta, thetadot')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()

# Plotting delta_v, delta_s, delta_m
plt.figure(figsize=(10, 6))

plt.plot(times_f, values_f, marker='.', linestyle=':', color='blue', linewidth=line_width, label='f')

plt.scatter(cp_times_f, cp_values_f, color='white', s=70, label='f Control Points', edgecolor='blue') 

plt.title('BeBOT Results: f')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()

# Show both plots at the end
plt.show()

