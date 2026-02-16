import pandas as pd
import numpy as np
import os

def analyze_species(filename, label):
    if not os.path.exists(filename):
        print(f"File {filename} not found.")
        return
    
    # Column names based on v_n_na_l v_n_cl_l v_n_na_r v_n_cl_r
    df = pd.read_csv(filename, sep='\s+', comment='#', 
                     names=['Step', 'Na_L', 'Cl_L', 'Na_R', 'Cl_R'])
    
    # Relative change from start
    dna_r = df['Na_R'].iloc[-1] - df['Na_R'].iloc[0]
    dcl_r = df['Cl_R'].iloc[-1] - df['Cl_R'].iloc[0]
    dna_l = df['Na_L'].iloc[-1] - df['Na_L'].iloc[0]
    dcl_l = df['Cl_L'].iloc[-1] - df['Cl_L'].iloc[0]
    
    # Net transfer (average of Gain in Right and Loss in Left)
    na_transfer = (dna_r - dna_l) / 2.0
    cl_transfer = (dcl_r - dcl_l) / 2.0
    
    print(f"\n=== Ion Selectivity Analysis: {label} ===")
    print(f"Na+ Net Transfer: {na_transfer:+.2f} ions")
    print(f"Cl- Net Transfer: {cl_transfer:+.2f} ions")
    
    charge_r = dna_r - dcl_r
    print(f"Induced Charge on Cold Side: {charge_r:+.2f}e")
    
    # Check for steady state by looking at the last 10% vs total
    recent_trend = df['Na_R'].iloc[-10:].mean() - df['Na_R'].iloc[0]
    print(f"Current Na+ Accumulation Trend: {recent_trend:+.2f}")

def analyze_all():
    print("--- Statistical Thermo-Ion Transport Analysis ---")
    
    cases = [
        ('funnel_window', 'Magic Window (Funnel)'),
        ('cylinder_control', 'Cylinder (Control)')
    ]
    
    for prefix, label in cases:
        all_dna_r = []
        all_dcl_r = []
        
        for i in range(4):
            filename = f'results/data.{prefix}_s{i}.dat'
            if not os.path.exists(filename): continue
            
            df = pd.read_csv(filename, sep='\s+', comment='#', 
                             names=['Step', 'Na_L', 'Cl_L', 'Na_R', 'Cl_R'])
            
            # Change over whole simulation
            dna_r = df['Na_R'].iloc[-1] - df['Na_R'].iloc[0]
            dcl_r = df['Cl_R'].iloc[-1] - df['Cl_R'].iloc[0]
            all_dna_r.append(dna_r)
            all_dcl_r.append(dcl_r)
            
        if not all_dna_r: continue
        
        na_mean = np.mean(all_dna_r)
        na_std = np.std(all_dna_r)
        cl_mean = np.mean(all_dcl_r)
        cl_std = np.std(all_dcl_r)
        charge_means = np.array(all_dna_r) - np.array(all_dcl_r)
        q_mean = np.mean(charge_means)
        q_std = np.std(charge_means)
        
        print(f"\n=== {label} Summary (N=4) ===")
        print(f"Na+ Transfer: {na_mean:+.2f} ± {na_std:.2f}")
        print(f"Cl- Transfer: {cl_mean:+.2f} ± {cl_std:.2f}")
        print(f"Net Charge Separation (Cold Side): {q_mean:+.2f} ± {q_std:.2f}e")

if __name__ == "__main__":
    analyze_all()
