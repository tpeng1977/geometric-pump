import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def load_multi_seed(prefix):
    files = glob.glob(f"results/data.{prefix}_s*.dat")
    all_data = []
    for f in files:
        try:
            df = pd.read_csv(f, sep='\s+', comment='#', 
                             names=['Step', 'Na_L', 'Cl_L', 'Na_R', 'Cl_R'])
            # Compute delta N (Right) - Initial N (Right) for each step
            df['dNa_R'] = df['Na_R'] - df['Na_R'].iloc[0]
            df['dCl_R'] = df['Cl_R'] - df['Cl_R'].iloc[0]
            df['dQ_R'] = df['dNa_R'] - df['dCl_R']
            all_data.append(df)
        except Exception as e:
            print(f"Skipping {f} due to error: {e}")
    return all_data

def plot_time_evolution():
    print("Generating time-evolution plot...")
    funnel_data = load_multi_seed('funnel_window')
    cylinder_data = load_multi_seed('cylinder_control')
    
    if not funnel_data or not cylinder_data:
        print("Missing data for time evolution.")
        return

    plt.figure(figsize=(10, 6))
    
    def get_stats(data_list):
        # We need to handle potentially different lengths if some crashed
        max_len = max(len(df) for df in data_list)
        # Re-index to the longest one and fill with NaN or just truncate
        # Here we truncate to the shortest to be safe for mean calculation
        min_len = min(len(df) for df in data_list)
        steps = data_list[0]['Step'].iloc[:min_len]
        q_matrix = np.array([df['dQ_R'].iloc[:min_len] for df in data_list])
        return steps, np.mean(q_matrix, axis=0), np.std(q_matrix, axis=0)

    f_steps, f_mean, f_std = get_stats(funnel_data)
    c_steps, c_mean, c_std = get_stats(cylinder_data)
    
    plt.plot(f_steps, f_mean, label='Magic Window (Funnel)', color='blue', linewidth=2)
    plt.fill_between(f_steps, f_mean - f_std, f_mean + f_std, color='blue', alpha=0.15)
    
    plt.plot(c_steps, c_mean, label='Control (Cylinder)', color='gray', linestyle='--', linewidth=2)
    plt.fill_between(c_steps, c_mean - c_std, c_mean + c_std, color='gray', alpha=0.15)
    
    plt.axhline(0, color='black', alpha=0.3)
    plt.title(f"Stochastic Charge Accumulation (N={len(funnel_data)} Seeds)", fontsize=14)
    plt.xlabel("Simulation Timestep", fontsize=12)
    plt.ylabel("Net Charge Separation $\Delta N$ (e)", fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    os.makedirs("figures", exist_ok=True)
    plt.savefig("figures/time_evolution_statistical.png", dpi=300, bbox_inches='tight')
    print("Saved: figures/time_evolution_statistical.png")

def plot_pdf_shift():
    print("Generating PDF shift histogram...")
    funnel_data = load_multi_seed('funnel_window')
    cylinder_data = load_multi_seed('cylinder_control')
    
    # FINAL WINDOW AVERAGING (Last 10,000 steps = 2 rows)
    f_finals = [df['dQ_R'].iloc[-2:].mean() for df in funnel_data if len(df) >= 2]
    c_finals = [df['dQ_R'].iloc[-2:].mean() for df in cylinder_data if len(df) >= 2]
    
    plt.figure(figsize=(9, 6))
    
    # Improved binning for N=20
    all_vals = f_finals + c_finals
    bins = np.linspace(min(all_vals) - 2, max(all_vals) + 2, 25)
    
    plt.hist(f_finals, bins=bins, alpha=0.6, label='Funnel (Experimental)', 
             color='blue', density=True, edgecolor='blue', linewidth=1.2)
    plt.hist(c_finals, bins=bins, alpha=0.3, label='Cylinder (Null Hypothesis)', 
             color='gray', density=True, edgecolor='black', linestyle='--', linewidth=1.0)
    
    plt.axvline(np.mean(f_finals), color='blue', linestyle='--', linewidth=2, label=f'Funnel Mean: {np.mean(f_finals):.2f}')
    plt.axvline(np.mean(c_finals), color='gray', linestyle='--', linewidth=2, label=f'Cylinder Mean: {np.mean(c_finals):.2f}')
    
    plt.title(f"Entropy PDF Reconstruction (N={len(funnel_data)} Seeds)", fontsize=14)
    plt.xlabel("Average Net Charge Shift $\Delta N$ (Last 10k Steps)", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.legend()
    plt.grid(True, axis='y', alpha=0.3)
    
    plt.savefig("figures/pdf_shift_reconstruction.png", dpi=300, bbox_inches='tight')
    print("Saved: figures/pdf_shift_reconstruction.png")

if __name__ == "__main__":
    plot_time_evolution()
    plot_pdf_shift()
