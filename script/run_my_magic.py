import os
import subprocess
from multiprocessing import Pool

LMP_EXE = "/opt/homebrew/bin/lmp_mpi"
SEEDS = [12345, 23456]

# Using templates based on the optimized in.funnel_run and in.cylinder_run
FUNNEL_TEMPLATE = """
# LAMMPS Input - Funnel Statistical Run
units           real
atom_style      full
boundary        p p p 
pair_style      lj/cut/coul/cut 10.0 10.0
dielectric      78.0

lattice         fcc 4.0
region          entire block -200 200 -200 200 -600 600 units box
create_box      2 entire

region          res_left    block -150 150 -150 150 -450  -40 units box
region          res_right   block -150 150 -150 150   40  450 units box
region          pore_hole   cone z 0 0 60 8 -60 60 units box
region          allowed     union 3 res_left res_right pore_hole

region          track_left  block -150 150 -150 150 -450 -70 units box
region          track_right block -150 150 -150 150   70  450 units box

region          spawn_left  block -140 140 -140 140 -400  -70 units box
region          spawn_right block -140 140 -140 140   70  400 units box

create_atoms    1 random 400 {seed1} spawn_left
create_atoms    2 random 400 {seed2} spawn_left
create_atoms    1 random 400 {seed3} spawn_right
create_atoms    2 random 400 {seed4} spawn_right

group           na_type      type 1
group           cl_type      type 2
group           hot_atoms    region res_left
group           cold_atoms   region res_right

mass            1 23.0
mass            2 35.5
set             type 1 charge 1.0
set             type 2 charge -1.0

pair_coeff      1 1 0.1 3.2
pair_coeff      1 2 0.1 3.0
pair_coeff      2 2 0.1 3.5

fix             w_all all wall/region allowed lj126 200.0 1.0 1.12
velocity        all create 300.0 {seed5}
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes
minimize        1.0e-4 1.0e-6 100 1000

fix             lang_left  hot_atoms  langevin 300.0 300.0 100.0 {seed6}
fix             lang_right cold_atoms langevin 300.0 300.0 100.0 {seed7}
fix             1 all nve

variable        n_na_l  equal count(na_type,track_left)
variable        n_cl_l  equal count(cl_type,track_left)
variable        n_na_r  equal count(na_type,track_right)
variable        n_cl_r  equal count(cl_type,track_right)

fix             ave_n all ave/time 1 1 5000 v_n_na_l v_n_cl_l v_n_na_r v_n_cl_r file results/{outfile}

timestep        0.25
run             4000000
"""

CYLINDER_TEMPLATE = FUNNEL_TEMPLATE.replace(
    "region          pore_hole   cone z 0 0 60 8 -60 60 units box",
    "region          pore_hole   cylinder z 0 0 34 -60 60 units box"
)

def run_single(args):
    name, template, i, seed, offset = args
    seed_final = seed + offset
    fname = f"in.{name}_s{i}"
    dname = f"data.{name}_s{i}.dat"
    lname = f"log.{name}_s{i}"
    
    content = template.format(
        seed1=seed_final, seed2=seed_final+1, seed3=seed_final+2, seed4=seed_final+3,
        seed5=seed_final+4, seed6=seed_final+5, seed7=seed_final+6,
        outfile=dname
    )
    
    os.makedirs("script", exist_ok=True)
    with open(f"script/{fname}", "w") as f:
        f.write(content)
        
    print(f"Starting {fname}...")
    subprocess.run(f"{LMP_EXE} -in script/{fname} > results/{lname} 2>&1", shell=True)
    print(f"Finished {fname}.")

def main():
    os.makedirs("results", exist_ok=True)
    
    tasks = []
    for i, seed in enumerate(SEEDS):
        tasks.append(("funnel_window", FUNNEL_TEMPLATE, i, seed, 0))
        tasks.append(("cylinder_control", CYLINDER_TEMPLATE, i, seed, 100))
    
    print(f"Planning 40 simulations (N=20 seeds) using parallel pool...")
    with Pool(processes=8) as pool:
        pool.map(run_single, tasks)

if __name__ == "__main__":
    main()
