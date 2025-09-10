mkdir -p results_leonardo/strong_scaling
mkdir -p results_leonardo/weak_scaling

sbatch Strong_scaling
sbatch Weak_scaling_1
sbatch Weak_scaling_2
sbatch Weak_scaling_4
sbatch Weak_scaling_8
