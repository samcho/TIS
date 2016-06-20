#ifndef PARAM_H
#define PARAM_H

void alloc_arrays();
void init_bonds(int);
void release_bonds();
void init_angles(int);
void release_angles();
void init_torsions(int);
void release_torsions();
void init_lj(int,int);
void release_lj();
void init_stacks(int);
void release_stacks();
void init_elec(int);
void release_elec();
void init_pos(int);
void release_pos();

void set_params(int);
void set_temp(double);

#endif
