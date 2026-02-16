import os
import subprocess
import time

# Configuration
LMP_EXE = "/opt/homebrew/bin/lmp_mpi"
VOLTAGES = [0.0, 0.001, 0.05] # Volts: 0 (No field), 1mV, 50mV (Reverse fields)
SEEDS = [12345, 23456, 34567, 45678]
STEPS = 200000
TIMESTEP = 0.1

TEMPLATE = """
# Power Sweep: Acetic Acid (Statistical Run)
# Voltage: {voltage} V, Seed: {seed}
units real
atom_style full
boundary p p f
pair_style lj/cut/coul/cut 12.0 12.0
dielectric 78.0

lattice fcc 4.0
region entire block -200 200 -200 200 -600 600 units box
create_box 2 entire
mass 1 1.008
mass 2 59.04

# Geometry
region l_res block -160 160 -160 160 -500 -45 units box
region r_res block -160 160 -160 160 45 500 units box
region pore cone z 0 0 100 4.5 -60 60 units box
region allowed union 3 l_res r_res pore
fix reflect all wall/reflect zlo -550 zhi 550

# Init
region spawn_l block -40 40 -40 40 -300 -200 units box
region spawn_r block -40 40 -40 40 200 300 units box
create_atoms 1 random 400 {seed} spawn_l
create_atoms 2 random 400 {seed_offset} spawn_l
create_atoms 1 random 400 {seed_offset_2} spawn_r
create_atoms 2 random 400 {seed_offset_3} spawn_r

group proton type 1
group acetate type 2
group mobile union proton acetate

set type 1 charge 1.0
set type 2 charge -1.0

pair_coeff 1 1 0.1 2.0
pair_coeff 1 2 0.1 4.5
pair_coeff 2 2 0.1 7.0

fix w_all mobile wall/region allowed lj126 1000.0 1.0 1.12

# Dynamics Init
velocity all create 350.0 {seed_offset_4}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
delete_atoms overlap 2.0 mobile all
minimize 1.0e-4 1.0e-6 1000 10000

# Apply Load Field
variable e_mag equal -1.0*{voltage}/1000.0
fix load_field all efield 0.0 0.0 v_e_mag

fix lang_mobile mobile langevin 350.0 350.0 100.0 {seed_offset_5}
fix 1 mobile nve

# Measurement
region track_l block -150 150 -150 150 -450 -100 units box
region track_r block -150 150 -150 150 100 450 units box

variable nhl equal count(proton,track_l)
variable nal equal count(acetate,track_l)
variable nhr equal count(proton,track_r)
variable nar equal count(acetate,track_r)

# Output
fix ave_n mobile ave/time 1 1 10000 v_nhl v_nal v_nhr v_nar file results/{outfile}

timestep {dt}
run {steps}
"""

def main():
    os.makedirs("results", exist_ok=True)
    os.makedirs("script", exist_ok=True)
    
    for v in VOLTAGES:
        for i, seed in enumerate(SEEDS):
            fname = f"in.power_sweep_{v}V_s{i}"
            dname = f"data.power_sweep_{v}V_s{i}.dat"
            lname = f"log.power_sweep_{v}V_s{i}"
            
            content = TEMPLATE.format(
                voltage=v, 
                seed=seed, 
                seed_offset=seed+1,
                seed_offset_2=seed+2,
                seed_offset_3=seed+3,
                seed_offset_4=seed+4,
                seed_offset_5=seed+123,
                outfile=dname, 
                dt=TIMESTEP, 
                steps=STEPS
            )
            
            with open(f"script/{fname}", "w") as f:
                f.write(content)
                
            print(f"Running {fname}...")
            subprocess.run(f"{LMP_EXE} -in script/{fname} > results/{lname} 2>&1", shell=True)
            time.sleep(0.5)

if __name__ == "__main__":
    main()
