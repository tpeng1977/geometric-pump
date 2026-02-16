import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import numpy as np

def analyze_sweep():
    print("--- Statistical Power Sweep Analysis ---")
    results = []
    
    # Analyze all data files in results/
    paths = glob.glob("results/data.power_sweep_*_s*.dat")
    
    for p in paths:
        try:
            # Filename format: data.power_sweep_0.0V_s0.dat
            parts = os.path.basename(p).split("_")
            voltage = float(parts[2].replace("V", ""))
            seed_idx = int(parts[3].replace("s", "").replace(".dat", ""))
            
            df = pd.read_csv(p, sep="\\s+", comment="#", names=["Step", "NHL", "NAL", "NHR", "NAR"])
            if df.empty: continue
            
            # Use last 50% for stability
            cutoff = len(df) // 2
            df_stable = df.iloc[cutoff:]
            
            avg_nhl = df_stable["NHL"].mean()
            avg_nhr = df_stable["NHR"].mean()
            delta_n = avg_nhr - avg_nhl
            
            results.append({
                "Voltage_V": voltage,
                "Seed": seed_idx,
                "DeltaN": delta_n
            })
        except Exception as e:
            print(f"Error processing {p}: {e}")
            
    if not results:
        print("No results found.")
        return

    res_df = pd.DataFrame(results)
    
    # Statistical Grouping
    summary = res_df.groupby("Voltage_V")["DeltaN"].agg(["mean", "std", "count"]).reset_index()
    summary["stderr"] = summary["std"] / np.sqrt(summary["count"])
    
    print("\n=== Power Sweep Statistical Summary (N=4 per voltage) ===")
    print(summary.to_string(index=False))
    
    summary.to_csv("results/power_sweep_statistical_analysis.csv", index=False)
    
    # Plotting
    plt.figure(figsize=(8, 5))
    plt.errorbar(summary["Voltage_V"] * 1000, summary["mean"], yerr=summary["stderr"], 
                 fmt='-o', capsize=5, color='darkgreen', linewidth=2, markersize=8, label='Net Charge Shift')
    
    plt.axhline(0, color='black', linestyle='--', alpha=0.3)
    plt.title("Statistical Work Extraction Profile (Acetic Acid)", fontsize=14)
    plt.xlabel("Retarding Potential (mV)", fontsize=12)
    plt.ylabel("$\Delta N$ (Steady State)", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    os.makedirs("figures", exist_ok=True)
    plt.savefig("figures/power_sweep_statistical.png", dpi=300)
    print("\nStatistical plot saved to: figures/power_sweep_statistical.png")

if __name__ == "__main__":
    analyze_sweep()
