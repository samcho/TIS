#include "global.h"
#include "random_generator.h"

int ncmd;
char cmd[mcmd+1][mwdsize];
char opt[mopt_tot+1][mwdsize];
int opt_ptr[mcmd+1];
char pathname[MAXPATHLEN];

// bonded info

double k_bnd; // bond spring constant
int nbnd; // number of bonds
int* ibead_bnd;
int* jbead_bnd;
double* pdb_dist; // pdb bond distances
int bnds_allocated = 0;
double e_bnd_coeff;

// angular info

double k_ang;
int nang;
int* ibead_ang;
int* jbead_ang;
int* kbead_ang;
double* pdb_ang;
int angs_allocated = 0;
double e_ang_coeff;

// torsion info

int ntor;
int* ibead_tor;
int* jbead_tor;
int* kbead_tor;
int* lbead_tor;
double* alpha;
double* beta;
double* gam;
double* delta;
double* cos_dih_ideal_dev;
double* sin_dih_ideal_dev;
int tors_allocated = 0;

// stacking info

int nstck;
rna_stck* rna_stcks;
int stcks_allocated = 0;
double alpha_st = 1.0;
double beta_st = 0.3;
double gamma_st = 1.0;

// rna-rna vdw

int ncon_att; // number of native contacts
int ncon_rep; // repulisve non-native contact
// neighbor list
int nnl_att;
int nnl_rep;

// pair list
int nil_att;
int nil_rep;

double coeff_att; // well-depth
double force_coeff_att;
double coeff_rep;
double sigma_rep;
double sigma_rep2;
double sigma_rep6;
double sigma_rep12;
double force_coeff_rep;
double rcut_nat;
int* ibead_lj_nat;
int* jbead_lj_nat;

double* lj_nat_pdb_dist;
double* lj_nat_pdb_dist2;
double* lj_nat_pdb_dist6;
double* lj_nat_pdb_dist12;

int* ibead_lj_non_nat;
int* jbead_lj_non_nat;

int* switch_fnb;

// neighbor / cell list
int* ibead_neighbor_list_att;
int* jbead_neighbor_list_att;

double* nl_lj_nat_pdb_dist;
double* nl_lj_nat_pdb_dist2;
double* nl_lj_nat_pdb_dist6;
double* nl_lj_nat_pdb_dist12;

int* ibead_neighbor_list_rep;
int* jbead_neighbor_list_rep;

// pair list
int* ibead_pair_list_att;
int* jbead_pair_list_att;

double* pl_lj_nat_pdb_dist;
double* pl_lj_nat_pdb_dist2;
double* pl_lj_nat_pdb_dist6;
double* pl_lj_nat_pdb_dist12;

int* ibead_pair_list_rep;
int* jbead_pair_list_rep;

int lj_rna_rna_allocated = 0;

// electrostatics

int nelec;
int* ibead_elec;
int* jbead_elec;
double ion_conc = 0.1; // in mols/L
double ion_str = ((ion_conc*z_p*z_p*Mol)+
		  (ion_conc*z_n*z_n*Mol))*elec*elec/2.0;;
double debye;
double felec_coeff;
double eelec_coeff;
int elec_allocated = 0;

// coordinates and associated params

int nbead;
coord* pos;
coord* unc_pos; // uncorrected positions
coord* vel;
coord* force;
coord* natpos;
int pos_allocated = 0;
int vel_allocated = 0;
int force_allocated = 0;
int natpos_allocated = 0;
int unc_pos_allocated = 0;

// native info

double nnc; // number of native contacts
int** ncmap;
int ncmap_allocated = 0;
double rcut_nat_plus_tol;
int* rna_base; // array which indicates whether or not a bead is a base
int rna_base_allocated;
int* rna_phosphate;
int rna_phosphate_allocated;

// miscellaneous run paramaters;

Ran_Gen generator; // random number generator
int run;
int restart = 0; // default is to start a new simulation
int rgen_restart = 0; // default don't restart random number generator
int sim_type = 1; // integration scheme; default is underdamped
double T; // temperature
int neighborlist = 0; // neighbor list cutoff method?
double minT; // min temperature determines crowder cutoffs
double boxl; // Length of an edge of the simulation box
double zeta; // friction coefficient
double nstep; // number of steps to take
double istep_restart = 0.0;
int nup;
int inlup;
int nnlup;
double h; // time step
double halfh;
double a1; // a1,a2,a3,a4,a5 are used for integration
double a2;
double a3;
double a4;
double a5;
char ufname[mwdsize+1];
char rcfname[mwdsize+1];
char cfname[mwdsize+1];
char unccfname[mwdsize+1];
char vfname[mwdsize+1];
char binfname[mwdsize+1];
char uncbinfname[mwdsize+1];
char iccnfigfname[mwdsize+1];
int binsave = 1; // default will save trajectory

// force and pot stuff

int nforce_term = 7; // ran,bnds,angs,tors,stacks,vdw,elec
int force_term_on[mforce_term+1] = { 0, 1, 1, 1, 1,
				     1, 1, 1, 0, 0, 0 };
force_term_Ptr force_term[mforce_term+1];

int npot_term = 6; // bnds,angs,tors,stacks,vdw,elec
int pot_term_on[mpot_term+1] = { 0, 1, 1, 1, 1,
				 1, 1, 0, 0, 0, 0 };
pot_term_Ptr pot_term[mpot_term+1];

//observables

double e_bnd,e_ang,e_tor,e_stack,e_elec;
double e_vdw_rr,e_vdw_rr_att,e_vdw_rr_rep;
double rna_etot,system_etot;
double chi;
double Q;
int contct_nat;
int contct_tot;
double end2endsq;
double rgsq;
double kinT;

