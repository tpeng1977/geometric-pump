import os
import subprocess
import numpy as np
import pandas as pd

# Configuration
LMP_EXE = "/opt/homebrew/bin/lmp_mpi"
N_SEEDS = 5
STEPS = 5000000
TIMESTEP = 0.1

TEMPLATES = {
    "isotope_h_d": """
# Isotope H/D Verification
units real
atom_style full
boundary p p f
pair_style lj/cut/coul/cut 10.0 10.0
dielectric 1.0
lattice fcc 4.0
region entire block -200 200 -200 200 -600 600 units box
create_box 2 entire
mass 1 1.008
mass 2 2.014
region l_res block -160 160 -160 160 -500 -45 units box
region r_res block -160 160 -160 160 45 500 units box
region pore cone z 0 0 100 4.0 -60 60 units box
region allowed union 3 l_res r_res pore
fix reflect all wall/reflect zlo -550 zhi 550
region spawn_l block -40 40 -40 40 -300 -200 units box
region spawn_r block -40 40 -40 40 200 300 units box
create_atoms 1 random 400 {seed1} spawn_l
create_atoms 2 random 400 {seed2} spawn_l
create_atoms 1 random 400 {seed3} spawn_r
create_atoms 2 random 400 {seed4} spawn_r
group h_type type 1
group d_type type 2
group mobile union h_type d_type
set type 1 charge 0.0
set type 2 charge 0.0
pair_coeff 1 1 0.1 3.0
pair_coeff 1 2 0.1 3.0
pair_coeff 2 2 0.1 3.0
fix w_all mobile wall/region allowed lj126 500.0 1.0 1.12
velocity all create 350.0 {seed5}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
minimize 1.0e-4 1.0e-6 100 1000
fix lang_all all langevin 350.0 350.0 100.0 {seed6}
fix 1 all nve
region track_l block -150 150 -150 150 -450 -100 units box
region track_r block -150 150 -150 150 100 450 units box
variable nhl equal count(h_type,track_l)
variable ndl equal count(d_type,track_l)
variable nhr equal count(h_type,track_r)
variable ndr equal count(d_type,track_r)
fix ave_n all ave/time 1 1 10000 v_nhl v_ndl v_nhr v_ndr file {output}
timestep {dt}
run {steps}
""",
    "acetic_low": """
# Acetic Acid Low Conc (n=400)
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
region l_res block -160 160 -160 160 -500 -45 units box
region r_res block -160 160 -160 160 45 500 units box
region pore cone z 0 0 100 4.5 -60 60 units box
region allowed union 3 l_res r_res pore
fix reflect all wall/reflect zlo -550 zhi 550
region spawn_l block -40 40 -40 40 -300 -200 units box
region spawn_r block -40 40 -40 40 200 300 units box
create_atoms 1 random 400 {seed1} spawn_l
create_atoms 2 random 400 {seed2} spawn_l
create_atoms 1 random 400 {seed3} spawn_r
create_atoms 2 random 400 {seed4} spawn_r
group proton type 1
group acetate type 2
group mobile union proton acetate
set type 1 charge 1.0
set type 2 charge -1.0
pair_coeff 1 1 0.1 2.0
pair_coeff 1 2 0.1 4.5
pair_coeff 2 2 0.1 7.0
fix w_all mobile wall/region allowed lj126 1000.0 1.0 1.12
velocity all create 350.0 {seed5}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
delete_atoms overlap 5.5 mobile all
minimize 1.0e-4 1.0e-6 1000 10000
fix lang_mobile mobile langevin 350.0 350.0 100.0 {seed6}
fix 1 mobile nve
region track_l block -150 150 -150 150 -450 -100 units box
region track_r block -150 150 -150 150 100 450 units box
variable nhl equal count(proton,track_l)
variable nal equal count(acetate,track_l)
variable nhr equal count(proton,track_r)
variable nar equal count(acetate,track_r)
fix ave_n mobile ave/time 1 1 10000 v_nhl v_nal v_nhr v_nar file {output}
timestep {dt}
run {steps}
""",
    "acetic_high": """
# Acetic Acid High Conc (n=1600)
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
region l_res block -160 160 -160 160 -500 -45 units box
region r_res block -160 160 -160 160 45 500 units box
region pore cone z 0 0 100 4.5 -60 60 units box
region allowed union 3 l_res r_res pore
fix reflect all wall/reflect zlo -550 zhi 550
region spawn_l block -40 40 -40 40 -300 -200 units box
region spawn_r block -40 40 -40 40 200 300 units box
create_atoms 1 random 1600 {seed1} spawn_l
create_atoms 2 random 1600 {seed2} spawn_l
create_atoms 1 random 1600 {seed3} spawn_r
create_atoms 2 random 1600 {seed4} spawn_r
group proton type 1
group acetate type 2
group mobile union proton acetate
set type 1 charge 1.0
set type 2 charge -1.0
pair_coeff 1 1 0.1 2.0
pair_coeff 1 2 0.1 4.5
pair_coeff 2 2 0.1 7.0
fix w_all mobile wall/region allowed lj126 1000.0 1.0 1.12
velocity all create 350.0 {seed5}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
delete_atoms overlap 5.5 mobile all
minimize 1.0e-4 1.0e-6 1000 10000
fix lang_mobile mobile langevin 350.0 350.0 100.0 {seed6}
fix 1 mobile nve
region track_l block -150 150 -150 150 -450 -100 units box
region track_r block -150 150 -150 150 100 450 units box
variable nhl equal count(proton,track_l)
variable nal equal count(acetate,track_l)
variable nhr equal count(proton,track_r)
variable nar equal count(acetate,track_r)
fix ave_n mobile ave/time 1 1 10000 v_nhl v_nal v_nhr v_nar file {output}
timestep {dt}
run {steps}
"""
}

def run_simulation(name, seeds, config_key):
    job_id = f"{config_key}_{name}"
    infile = f"script/in.verify_{job_id}"
    outfile = f"results/data.verify_{job_id}.dat"
    logfile = f"results/log.verify_{job_id}"

    seed_dict = {f"seed{i+1}": s for i, s in enumerate(seeds)}

    content = TEMPLATES[config_key].format(
        steps=STEPS, dt=TIMESTEP, output=outfile, **seed_dict
    )

    os.makedirs("script", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    with open(infile, "w") as f:
        f.write(content)

    print(f"Running {job_id}...")
    cmd = f"mpirun -np 3 {LMP_EXE} -in {infile} > {logfile} 2>&1"
    subprocess.run(cmd, shell=True, cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return outfile

def analyze_isotope(files):
    results = []
    for f in files:
        if not os.path.exists(f):
            continue
        try:
            df = pd.read_csv(f, sep=r"\s+", comment="#", names=["Step", "H_L", "D_L", "H_R", "D_R"])
            if df.empty:
                continue
            last = df.iloc[-1]
            delta = (last["D_R"] - last["D_L"]) - (last["H_R"] - last["H_L"])
            results.append(delta)
        except Exception as e:
            print(f"Error reading {f}: {e}")
    return results

def analyze_acetic(files):
    results = []
    for f in files:
        if not os.path.exists(f):
            continue
        try:
            df = pd.read_csv(f, sep=r"\s+", comment="#", names=["Step", "H_L", "OAc_L", "H_R", "OAc_R"])
            if df.empty:
                continue
            last = df.iloc[-1]
            delta = last["H_R"] - last["H_L"]
            results.append(delta)
        except Exception as e:
            print(f"Error reading {f}: {e}")
    return results

def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)

    report = "# Isothermal Rectification Verification Report\n\n"

    seed_sets = []
    for i in range(N_SEEDS):
        seed_sets.append(np.random.randint(1000, 99999, size=6))

    # 1. Isotope
    files = []
    for i, seeds in enumerate(seed_sets):
        files.append(run_simulation(f"seed{i}", seeds, "isotope_h_d"))
    res = analyze_isotope(files)
    mean = np.mean(res) if res else float("nan")
    sem = np.std(res) / np.sqrt(len(res)) if len(res) > 1 else 0.0
    report += "## Isotope H/D (Momentum Filter)\n"
    report += f"- N={len(res)}, Steps={STEPS}\n"
    report += f"- Mean Separation (D_tip - H_tip): **{mean:.2f} +/- {sem:.2f}**\n\n"

    # 2. Acetic Low
    files = []
    for i, seeds in enumerate(seed_sets):
        files.append(run_simulation(f"seed{i}", seeds, "acetic_low"))
    res = analyze_acetic(files)
    mean = np.mean(res) if res else float("nan")
    sem = np.std(res) / np.sqrt(len(res)) if len(res) > 1 else 0.0
    report += "## Acetic Acid Low Conc (Entropic Pumping)\n"
    report += f"- N={len(res)}, Steps={STEPS}\n"
    report += f"- Mean H+ Shift (Right - Left): **{mean:.2f} +/- {sem:.2f}**\n\n"

    # 3. Acetic High
    files = []
    for i, seeds in enumerate(seed_sets):
        files.append(run_simulation(f"seed{i}", seeds, "acetic_high"))
    res = analyze_acetic(files)
    mean = np.mean(res) if res else float("nan")
    sem = np.std(res) / np.sqrt(len(res)) if len(res) > 1 else 0.0
    report += "## Acetic Acid High Conc (Jamming Switch)\n"
    report += f"- N={len(res)}, Steps={STEPS}\n"
    report += f"- Mean H+ Shift (Right - Left): **{mean:.2f} +/- {sem:.2f}**\n\n"

    with open("results/verification_report.md", "w") as f:
        f.write(report)
    print("Verification complete. Report written to results/verification_report.md")

if __name__ == "__main__":
    main()
