import os
import pandas as pd

def update_or_create_csv():
    # File paths
    output_file = "data_general.csv"
    theta_file = "theta.csv"
    thetadot_file = "thetadot.csv"
    f_file = "f.csv"
    
    # Columns for the output CSV
    columns = ["m_value", "M_value", "L_value", "g_value", "d_value", "b_value", 
               "theta_cmd", "theta_0", "thetadot_cmd", "theta_dot_0", "F_value", "t0"]
    
    # Read the required data from the source files
    theta_df = pd.read_csv(theta_file)
    thetadot_df = pd.read_csv(thetadot_file)
    f_df = pd.read_csv(f_file)

    # Extract required data starting from the first line after the title
    t0_values = theta_df["Time"].tolist()
    theta_0_values = theta_df["Value"].tolist()
    theta_dot_0_values = thetadot_df["Value"].tolist()
    f_values = f_df["Value"].tolist()
    
    # Determine the length of data
    num_rows = len(t0_values)
    
    # Prepare default data for other columns
    default_zeros = [0] * num_rows
    
    # Create a DataFrame with the desired structure
    data = {
        "m_value": default_zeros,
        "M_value": default_zeros,
        "L_value": default_zeros,
        "g_value": default_zeros,
        "d_value": default_zeros,
        "b_value": default_zeros,
        "theta_cmd": default_zeros,
        "theta_0": theta_0_values,
        "thetadot_cmd": default_zeros,
        "theta_dot_0": theta_dot_0_values,
        "F_value": f_values,
        "t0": t0_values,
    }
    new_data_df = pd.DataFrame(data)

    # Check if the output file exists
    if not os.path.exists(output_file):
        # File does not exist; create it and write the header
        new_data_df.to_csv(output_file, index=False, sep="\t")
        print(f"{output_file} created and filled with data.")
    else:
        # File exists; append data
        with open(output_file, "r") as f:
            lines = f.readlines()

        # Find the last filled line
        filled_line_index = len(lines) - 1
        while filled_line_index > 0 and lines[filled_line_index].strip() == "":
            filled_line_index -= 1
        
        # Insert the title line below the last filled line
        new_lines = lines[:filled_line_index + 1]
        new_lines.append("\t".join(columns) + "\n")
        
        # Write back the file with the new header
        with open(output_file, "w") as f:
            f.writelines(new_lines)
        
        # Append the new data below the new title line
        new_data_df.to_csv(output_file, mode='a', index=False, header=False, sep="\t")
        print(f"Data appended to {output_file}.")

# Run the function
update_or_create_csv()
