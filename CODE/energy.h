#ifndef ENERGY_H
#define ENERGY_H

void set_potential();
void set_forces();
void clear_forces();
void energy_eval();
void force_eval();

void random_force();
void bond_energy();
void bond_forces();
void angular_energy();
void angular_forces();
void torsional_energy();
void torsional_forces();

void vdw_energy();
void vdw_forces();
void stacking_energy();
void stacking_forces();
void electrostatic_energy();
void electrostatic_forces();

#endif /* ENERGY_H */
