import pandas as pd
import matplotlib.pyplot as plt

# File path to the CSV file
file_path = "data_general.csv"

# Initialize lists for plotting
theta_0_values = []
theta_dot_0_values = []
row_indices = []

try:
    # Load the file with automatic delimiter detection
    data = None
    # Load the file with the correct delimiter
    try:
        data = pd.read_csv(file_path, delimiter=',')  # Explicitly specify comma-delimited
        print("Loaded file as comma-delimited.")
    except Exception as e:
        print(f"Error loading file: {e}")


    # Verify the loaded data
    print("Data preview:")
    print(data.head())

    # Check if 't0', 'theta_0', and 'theta_dot_0' columns exist
    if 't0' not in data.columns or 'theta_0' not in data.columns or 'theta_dot_0' not in data.columns:
        print("Required columns ('t0', 'theta_0', 'theta_dot_0') are missing from the file.")
    else:
        # Loop through rows to find where 't0' is 0 and extract thetas
        for i, row in data.iterrows():
            try:
                if float(row['t0']) == 0:  # Check if 't0' is 0
                    theta_0_values.append(float(row['theta_0']))
                    theta_dot_0_values.append(float(row['theta_dot_0']))
                    row_indices.append(i)  # Store row index for debugging
                    print(f"Row {i}: Found theta_0={row['theta_0']}, theta_dot_0={row['theta_dot_0']}, t0={row['t0']}")
            except (ValueError, KeyError):
                print(f"Skipping invalid or incomplete data at row {i}")

    # Check if any data was extracted
    if len(theta_0_values) == 0 or len(theta_dot_0_values) == 0:
        print("No valid data found for plotting.")
    else:
        # Create a scatter plot
        plt.figure(figsize=(10, 6))
        scatter = plt.scatter(
            theta_0_values,
            theta_dot_0_values,
            c=row_indices,
            cmap='viridis',
            edgecolor='k',
            s=100
        )
        plt.colorbar(scatter, label='Row Index')
        plt.xlabel('theta_0')
        plt.ylabel('theta_dot_0')
        plt.title('Theta_0 vs Theta_Dot_0 (Rows with t0=0)')
        plt.grid(True)
        plt.show()

except FileNotFoundError:
    print(f"File not found: {file_path}")
except pd.errors.ParserError as e:
    print(f"Error parsing the file: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
