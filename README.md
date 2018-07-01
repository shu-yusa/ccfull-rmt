# Coupled-channels calculation with noncollective excitations
Code for calculating quasi-elastic scattering and fusion reaction
using the coupled-channels method.
In addition to collective excitations such as vibrations and rotations, this takes into accout
noncollective excitations using random matrix theory.

## Prerequisites
`make` and `gfortran` need to be installed. `mpich` enables distributed calculation.
```bash
brew install gfortran
brew install mpich
```

## Compile & Run
Compile the source codes by
```bash
make
```

This will produce an executable file `a.out`.  
Run calculation with a single node by
```bash
./a.out
```

To run on multiple nodes, specify nodes in `machinefile` and run
```bash
./exe.sh `number of nodes` a.out
```

## Input file
* `input_scat_exact_noncoll` ... Input file for the calculation. You can configure reaction details and output directory.
* `input_rmt` ... Parameters for random matrix.

## Details of source codes and modules
* `global_constant.f90` ... Definition of constants.
* `input_data.f90` ... Read input parameters from the input file.
* `scat_noncoll.f90` ... The main program.
* `potentials.f90` ... Potential class is defined. Methods for calculating Coulomb, nuclear, and centrifugal potentials are defined.
* `relative_potential.f90` ... Class for calculating the relative potential for the reaction.
* `coupling_matrix.f90` ... Class for calculating the coupling matrix for the reaction.
* `coupled_channels.f90` ... Class for calculating coupled-channels calculation by Numerov method. 
* `calc_profile.f90` ... Class for output calculation details into a file `calc_info`.
* `utils/` ... Utility functions are contained in this directory.

## Output
Calculation results are output in the directory specified in the input file. 
Outputs are consisted of the following files.
* `calc_info` ... Calculation configurations.
* `fusion.dat` ... Fusion cross sections.
* `fus_bar_dist.dat` ... Fusion barrier distributions.
* `pot_nucl.dat` ... Nuclear potentions.
* `potential.dat` ... Real and imaginary parts of the total potential.
* `angular/` ... Quasi-elastic scattering cross sections and barrier distributions at specified angles.
* `angular_dist/` ... Angular distribution of elastic differential scattering cross sections.
* `Q_val_dist/` ... Q-value distribution of the scattering.

## References
[Tohoku University Repository（TOUR）](https://tohoku.repo.nii.ac.jp/?action=pages_view_main&active_action=repository_view_main_item_detail&item_id=70405&item_no=1&page_id=33&block_id=38)
