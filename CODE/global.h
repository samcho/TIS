#ifndef GLOBAL_H
#define GLOBAL_H

#include <cmath>
#include <cstring>

const int mcmd = 100; // maximum number of input commands
const int mopt = 10; // maximum number of options associated with a command
const int mopt_tot = mcmd*mopt; // max total number of options
const int mwdsize = 1024; // maximum number of characters in a word
const size_t MAXPATHLEN = 2048;

extern char cmd[][mwdsize];
extern char opt[][mwdsize]; // holds the options
extern int opt_ptr[]; // holds index of first option correspond to given cmd
extern int ncmd; // number of input commands
extern int nopt_tot; // total # of options
extern char pathname[];

// bonded info

extern double k_bnd;
extern int nbnd;
extern int* ibead_bnd;
extern int* jbead_bnd;
extern double* pdb_dist;
extern int bnds_allocated;
extern double e_bnd;
extern double e_bnd_coeff;

// angular info

extern double k_ang;
extern int nang;
extern int* ibead_ang;
extern int* jbead_ang;
extern int* kbead_ang;
extern double* pdb_ang;
extern int angs_allocated;
extern double e_ang;
extern double e_ang_coeff;

// torsion info

extern int ntor;
extern int* ibead_tor;
extern int* jbead_tor;
extern int* kbead_tor;
extern int* lbead_tor;
extern double* alpha;
extern double* beta;
extern double* gam;
extern double* delta;
extern double* cos_dih_ideal_dev;
extern double* sin_dih_ideal_dev;
extern int tors_allocated;
extern double e_tor;

// stacking info

const int nstck_ang = 4;
const int nstck_dist = 2;
const int nstck_tor = 2;

class rna_stck {

 public:
  rna_stck() { };
  ~rna_stck() { };
  int ibead_ang[nstck_ang+1];
  int jbead_ang[nstck_ang+1];
  int kbead_ang[nstck_ang+1];
  double pdb_ang[nstck_ang+1];
  double cos_pdb_ang[nstck_ang+1];
  double sin_pdb_ang[nstck_ang+1];
  int ibead_dist[nstck_dist+1];
  int jbead_dist[nstck_dist+1];
  double pdb_dist[nstck_dist+1];
  int ibead_tor[nstck_tor+1];
  int jbead_tor[nstck_tor+1];
  int kbead_tor[nstck_tor+1];
  int lbead_tor[nstck_tor+1];
  double pdb_tor[nstck_tor+1];
  double cos_pdb_tor[nstck_tor+1];
  double sin_pdb_tor[nstck_tor+1];
  double delta_H;
  double delta_S;
  double delta_G; // make sure to change this if change temp

};  

extern int nstck;
extern rna_stck* rna_stcks;
extern int stcks_allocated;
extern double alpha_st;
extern double beta_st;
extern double gamma_st;
extern double e_stack;

// rna-rna vdw info

extern int ncon_att; // number of native contacts
extern int ncon_rep; // repulisve non-native contact

// neighbor list
extern int nnl_att;
extern int nnl_rep;

// pair list
extern int nil_att;
extern int nil_rep;

extern double coeff_att; // well-depth
extern double force_coeff_att;
extern double coeff_rep;
extern double force_coeff_rep;
extern double sigma_rep;
extern double sigma_rep2;
extern double sigma_rep6;
extern double sigma_rep12;
extern double rcut_nat; // dist cutoff for defining nbcontact via potential
extern int* ibead_lj_nat;
extern int* jbead_lj_nat;
extern double* lj_nat_pdb_dist;
extern double* lj_nat_pdb_dist2; // 2nd power of the pdb distance
extern double* lj_nat_pdb_dist6; // 6th power of the pdb distance
extern double* lj_nat_pdb_dist12; // 12th power of the pdb distance
extern int* ibead_lj_non_nat;
extern int* jbead_lj_non_nat;
extern int lj_rna_rna_allocated;
extern int* switch_fnb;
extern double e_vdw_rr;
extern double e_vdw_rr_att;
extern double e_vdw_rr_rep;

// electrostatic info

extern int nelec;
extern int* ibead_elec;
extern int* jbead_elec;
extern int elec_allocated;
extern double debye;
const double NA = 6.02e23;
const double Jtocal = 0.23923;
const double kilo = 0.001; // i.e. kJ/J, kcal/cal, etc.
const double z_e = 1.0; // negative charge on phosphate group
// const double epsilon_r = 10.0; // dielectric constant
const double epsilon_r = 5.0; // dielectric constant
const double kelec = 8.99e9; // 1/(4*pi*epsilon_0) in JmC^{-2}
const double elec = 1.602e-19; // electron charge
const double z_p = 1.0; // Na+
const double z_n = -1.0; // Cl-
// const double ion_conc = 0.2; // in mols/L
//extern double ion_conc = 0.1; // in mols/L
extern double ion_conc; // in mols/L
const double Mol = 6.02e23/0.001;
extern double ion_str;
const double meter2angstrom = 1.0e10; // angstroms/meter
const double Debye_limit = 0.5;
const double JtoKCpm = NA*Jtocal*kilo;
const double KCpmtoJ = 6.9412e-21;
extern double felec_coeff;
extern double eelec_coeff;
extern double e_elec;

// neighbor / cell list
extern int* ibead_neighbor_list_att;
extern int* jbead_neighbor_list_att;

extern double* nl_lj_nat_pdb_dist;
extern double* nl_lj_nat_pdb_dist2;
extern double* nl_lj_nat_pdb_dist6;
extern double* nl_lj_nat_pdb_dist12;

extern int* ibead_neighbor_list_rep;
extern int* jbead_neighbor_list_rep;

// pair list
extern int* ibead_pair_list_att;
extern int* jbead_pair_list_att;

extern double* pl_lj_nat_pdb_dist;
extern double* pl_lj_nat_pdb_dist2;
extern double* pl_lj_nat_pdb_dist6;
extern double* pl_lj_nat_pdb_dist12;

extern int* ibead_pair_list_rep;
extern int* jbead_pair_list_rep;

extern int lj_rna_rna_allocated;
extern int* switch_fnb;
extern double e_vdw_rr;
extern double e_vdw_rr_att;
extern double e_vdw_rr_rep;

// coordinates and associated params

class coord {

 public:
  coord() {
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }
  ~coord() { };
  double x;
  double y;
  double z;

};

extern int nbead;
extern int ncrowder;
extern int nbead_tot;
extern coord* pos;
extern coord* unc_pos;
extern coord* vel;
extern coord* force;
extern coord* natpos; // native position vectors
extern int pos_allocated;
extern int unc_pos_allocated;
extern int vel_allocated;
extern int force_allocated;
extern int natpos_allocated;

// native info

extern double nnc; // number of native contacts
extern int** ncmap;
extern int ncmap_allocated;
extern double rcut_nat_plus_tol;
extern double** native_distance;
extern int native_distance_allocated;
extern int* rna_base;
extern int rna_base_allocated;
extern int* rna_phosphate;
extern int rna_phosphate_allocated;

// miscellaneous run paramaters;

extern int run;
extern class Ran_Gen generator; // the random number generator
extern int restart; // are we restarting an old simulation?
extern int rgen_restart; // should we restart the random number generator?
extern int sim_type; // integration scheme 1 = underdamped; 2 = overdamped
extern double T; // temperature (kcal/mol)
extern int neighborlist; // neighbor list cutoff method?
extern double minT; // minimum temperature determines crowder cutoffs
extern double boxl; // Length of an edge of the simulation box
extern double zeta; // friction coefficient
extern double nstep; // number of steps to take
extern double istep_restart; // which step to we restart from?
extern int nup;
extern int inlup;
extern int nnlup;
extern double h; // time step
extern double halfh;
extern double a1; // a1,a2,a3,a4 are used for integration
extern double a2;
extern double a3;
extern double a4;
extern double a5;
extern char ufname[];
extern char rcfname[];
extern char cfname[];
extern char unccfname[];
extern char vfname[];
extern char binfname[];
extern char uncbinfname[];
extern char iccnfigfname[];
extern int binsave;

const double pi = acos(-1);

/* potential stuff */

extern int npot_term;// number of terms in the potential
const int mpot_term = 10;// max number of terms in potential
extern int pot_term_on[];// is a particular term in the potential on?
typedef void (*pot_term_Ptr) ();
/* array of pointers to functions;
   each element is for evaluating a 
   particular term in the potential */
extern pot_term_Ptr pot_term[]; 
extern double rna_etot;
extern double system_etot;

/* force stuff */

extern int nforce_term; // number of force terms
const int mforce_term = 10; // max number of force terms
extern int force_term_on[]; // is a particular force term on?
typedef void (*force_term_Ptr) ();
extern force_term_Ptr force_term[]; // array of pointers to functions -- each elements is for evaluating a particular type of force


// observables

extern double chi;
extern double Q;
extern int contct_nat;
extern int contct_tot;
extern double end2endsq;
extern double rgsq;
extern double kinT;

// conversion factors;

const double kcalpmol2K = 503.15;

#endif /* GLOBAL_H */
