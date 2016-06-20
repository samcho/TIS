#ifndef TIS_H
#define TIS_H

#include "global.h"

void ex_cmds(); // sequentially execute cmds found in input_file
void simulation_ctrl();
void underdamped_ctrl();
void overdamped_ctrl();
void underdamped_iteration(coord*);
void overdamped_iteration(coord*);
void calculate_observables(coord*);
void print_sim_params();

void update_neighbor_list();
void update_pair_list();

#endif /* TIS_H */
