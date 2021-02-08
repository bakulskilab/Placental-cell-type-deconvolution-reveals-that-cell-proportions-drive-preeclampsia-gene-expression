# Large Memory Slurm Interactive Job Submission
# Standard Memeory Nodes go up to 7GB
# LargeMem goes up to 41.75

# Identify running jobs under my username
squeue -u kyleac

# getting weird errors with this setup
# https://groups.google.com/forum/#!topic/slurm-users/nco5m8SIgE0 suggests maybe ntasks the issue?
srun --nodes=6 --account=bakulski1 --ntasks-per-node=1 --mem-per-cpu=6GB --pty /bin/bash

# re-run with ntasks = 4, default shown on https://arc-ts.umich.edu/greatlakes/slurm-user-guide/
srun --nodes=6 --account=bakulski1 --ntasks-per-node=4 --mem-per-cpu=6GB --pty /bin/bash