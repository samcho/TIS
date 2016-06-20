#include <iostream>
#include <cstring>
#include "param.h"
#include "global.h"
#include "random_generator.h"

void alloc_arrays()
{

  using namespace std;

  // bonds

  k_bnd = 20.0;
  e_bnd_coeff = k_bnd/2.0;
  nbnd = 226;
  ibead_bnd = new int[nbnd+1];
  jbead_bnd = new int[nbnd+1];
  pdb_dist = new double[nbnd+1];
  bnds_allocated = 1;

  // angles

  k_ang = 20.0;
  e_ang_coeff = k_ang/2.0;
  nang = 149;
  ibead_ang = new int[nang+1];
  jbead_ang = new int[nang+1];
  kbead_ang = new int[nang+1];
  pdb_ang = new double[nang+1];
  angs_allocated = 1;

  // torsions

  ntor = 148;
  ibead_tor = new int[ntor+1];
  jbead_tor = new int[ntor+1];
  kbead_tor = new int[ntor+1];
  lbead_tor = new int[ntor+1];
  alpha = new double[ntor+1];
  beta = new double[ntor+1];
  gam = new double[ntor+1];
  delta = new double[ntor+1];
  cos_dih_ideal_dev = new double[ntor+1];
  sin_dih_ideal_dev = new double[ntor+1];
  tors_allocated = 1;

  // stacks

  nstck = 17;
  rna_stcks = new rna_stck[nstck+1];
  stcks_allocated = 1;

  // rna-rna vdw

  ncon_att = 147;
  ncon_rep = 25129;
  // neighbor list
  nnl_att = 0;
  nnl_rep = 0;
  // pair list
  nil_att = 0;
  nil_rep = 0;

  rcut_nat = 7.0;
  coeff_att = 1.5;
  force_coeff_att = -12.0*coeff_att;
  coeff_rep = 1.0;
  sigma_rep = 3.4;
  sigma_rep2 = sigma_rep*sigma_rep;
  sigma_rep6 = sigma_rep2*sigma_rep2*sigma_rep2;
  sigma_rep12 = sigma_rep6*sigma_rep6;
  force_coeff_rep = -6.0*coeff_rep;
  ibead_lj_nat = new int[ncon_att+1];
  jbead_lj_nat = new int[ncon_att+1];
  lj_nat_pdb_dist = new double[ncon_att+1];
  lj_nat_pdb_dist2 = new double[ncon_att+1];
  lj_nat_pdb_dist6 = new double[ncon_att+1];
  lj_nat_pdb_dist12 = new double[ncon_att+1];
  ibead_lj_non_nat = new int[ncon_rep+1];
  jbead_lj_non_nat = new int[ncon_rep+1];
  switch_fnb = new int[ncon_rep+1];
  lj_rna_rna_allocated = 1;

  // electrostatics

  nelec = 2775;
  felec_coeff = -z_e*z_e*elec*elec*kelec/epsilon_r*JtoKCpm*meter2angstrom;
  eelec_coeff = z_e*z_e*elec*elec*kelec/epsilon_r*JtoKCpm*meter2angstrom;
  ibead_elec = new int[nelec+1];
  jbead_elec = new int[nelec+1];
  elec_allocated = 1;

  // coordinates
  
  nbead = 227;
  rcut_nat_plus_tol = 7.5;
  pos = new coord[nbead+1];
  unc_pos = new coord[nbead+1];
  vel = new coord[nbead+1];
  force = new coord[nbead+1];
  rna_base = new int [nbead+1];
  rna_phosphate = new int [nbead+1];
  pos_allocated = 1;
  unc_pos_allocated = 1;
  vel_allocated = 1;
  force_allocated = 1;
  rna_base_allocated = 1;
  rna_phosphate_allocated = 1;

  // miscellaneous run parameters

  run = 1;
  generator.set_seed(-100-run);
  T = 0.6; // kcal/mol
  minT = 0.38;
  boxl = 200.0;
  zeta = 5.0e-2; // 0.05*tau^{-1} = friction coeff
  nstep = 5e7;
  nup = 1000;
  h = 2.5e-3;
  halfh = h/2.0;
  a1 = h*(1.0-zeta*halfh);
  a2 = h*halfh;
  a3 = (1.0-h*zeta/2.0+(h*zeta)*(h*zeta)/4.0)/h;
  a4 = halfh*(1.0-h*zeta/2.0);
  a5 = h/zeta;
  strcpy(ufname,"update.out");
  strcpy(rcfname,"restart_c.dat");
  strcpy(cfname,"coord.out");
  strcpy(unccfname,"unccoord.out");
  strcpy(vfname,"veloc.out");
  strcpy(binfname,"traj.bin");
  strcpy(uncbinfname,"traj_uncorrected.bin");
  strcpy(iccnfigfname,"iccnfig.xyz");
  debye = sqrt(epsilon_r*T*KCpmtoJ/(8.0*pi*kelec*ion_str))*meter2angstrom;
  if( debye < Debye_limit ) debye = Debye_limit;

}

void init_bonds(int numbonds)
{
 
  using namespace std;

  nbnd = numbonds;
  ibead_bnd = new int[numbonds+1];
  jbead_bnd = new int[numbonds+1];
  pdb_dist = new double[numbonds+1];
  bnds_allocated = 1;

}

void release_bonds()
{

  using namespace std;

  delete [] ibead_bnd;
  delete [] jbead_bnd;
  delete [] pdb_dist;
  bnds_allocated = 0;

}

void init_angles(int numangs)
{
 
  using namespace std;

  nang = numangs;
  ibead_ang = new int[numangs+1];
  jbead_ang = new int[numangs+1];
  kbead_ang = new int[numangs+1];
  pdb_ang = new double[numangs+1];
  angs_allocated = 1;

}

void release_angles()
{

  using namespace std;

  delete [] ibead_ang;
  delete [] jbead_ang;
  delete [] kbead_ang;
  delete [] pdb_ang;
  angs_allocated = 0;

}

void init_torsions(int numtors)
{
 
  using namespace std;

  ntor = numtors;
  ibead_tor = new int[numtors+1];
  jbead_tor = new int[numtors+1];
  kbead_tor = new int[numtors+1];
  lbead_tor = new int[numtors+1];
  alpha = new double[numtors+1];
  beta = new double[numtors+1];
  gam = new double[numtors+1];
  delta = new double[numtors+1];
  cos_dih_ideal_dev = new double[numtors+1];
  sin_dih_ideal_dev = new double[numtors+1];
  tors_allocated = 1;

}

void release_torsions()
{

  using namespace std;

  delete [] ibead_tor;
  delete [] jbead_tor;
  delete [] kbead_tor;
  delete [] lbead_tor;
  delete [] alpha;
  delete [] beta;
  delete [] gam;
  delete [] delta;
  delete [] cos_dih_ideal_dev;
  delete [] sin_dih_ideal_dev;
  tors_allocated = 0;

}

void init_lj(int numatt, int numrep ) 
{

  using namespace std;

  ncon_att = numatt;
  ncon_rep = numrep;
  ibead_lj_nat = new int[numatt+1];
  jbead_lj_nat = new int[numatt+1];
  lj_nat_pdb_dist = new double[numatt+1];
  lj_nat_pdb_dist2 = new double[numatt+1];
  lj_nat_pdb_dist6 = new double[numatt+1];
  lj_nat_pdb_dist12 = new double[numatt+1];
  ibead_lj_non_nat = new int[numrep+1];
  jbead_lj_non_nat = new int[numrep+1];
  switch_fnb = new int[numrep+1];

  ibead_neighbor_list_att = new int[numatt+1];
  jbead_neighbor_list_att = new int[numatt+1];
  nl_lj_nat_pdb_dist = new double[numatt+1];
  nl_lj_nat_pdb_dist2 = new double[numatt+1];
  nl_lj_nat_pdb_dist6 = new double[numatt+1];
  nl_lj_nat_pdb_dist12 = new double[numatt+1];
  ibead_neighbor_list_rep = new int[numrep+1];
  jbead_neighbor_list_rep = new int[numrep+1];

  ibead_pair_list_att = new int[numatt+1];
  jbead_pair_list_att = new int[numatt+1];
  pl_lj_nat_pdb_dist = new double[numatt+1];
  pl_lj_nat_pdb_dist2 = new double[numatt+1];
  pl_lj_nat_pdb_dist6 = new double[numatt+1];
  pl_lj_nat_pdb_dist12 = new double[numatt+1];
  ibead_pair_list_rep = new int[numrep+1];
  jbead_pair_list_rep = new int[numrep+1];

  lj_rna_rna_allocated = 1;

}

void release_lj()
{

  using namespace std;

  delete [] ibead_lj_nat;
  delete [] jbead_lj_nat;
  delete [] lj_nat_pdb_dist;
  delete [] lj_nat_pdb_dist2;
  delete [] lj_nat_pdb_dist6;
  delete [] lj_nat_pdb_dist12;
  delete [] ibead_lj_non_nat;
  delete [] jbead_lj_non_nat;
  delete [] switch_fnb;
  lj_rna_rna_allocated = 0;

}

void init_stacks(int numstcks)
{

  using namespace std;

  nstck = numstcks;
  rna_stcks = new rna_stck[numstcks+1];
  stcks_allocated = 1;

}

void release_stacks()
{

  using namespace std;

  delete [] rna_stcks;
  stcks_allocated = 0;

}

void init_elec(int numelec)
{

  using namespace std;

  nelec = numelec;
  ibead_elec = new int[numelec+1];
  jbead_elec = new int[numelec+1];
  elec_allocated = 1;

}

void release_elec()
{
  
  using namespace std;

  delete [] ibead_elec;
  delete [] jbead_elec;
  elec_allocated = 0;
  
}

void init_pos(int nbead)
{

  using namespace std;

  unc_pos = new coord[nbead+1];
  pos = new coord[nbead+1];

  vel = new coord[nbead+1];
  force = new coord[nbead+1];

  pos_allocated = 1;
  unc_pos_allocated = 1;
  vel_allocated = 1;
  force_allocated = 1;
}

void release_pos()
{

  using namespace std;

  delete [] unc_pos;
  delete [] pos;

  delete [] vel;
  delete [] force;

  pos_allocated = 0;
  unc_pos_allocated = 0;
  vel_allocated = 0;
  force_allocated = 0;
}

void set_params(int icmd)
{

  using namespace std;
  char oline[1024];
  int iopt;

  if( !strcmp(opt[opt_ptr[icmd]],"dynamics") ) { // set the type of simulation
    if( !strcmp(opt[opt_ptr[icmd]+1],"underdamped") ) {
      sim_type = 1; // low-friction limit for sampling
    } else if( !strcmp(opt[opt_ptr[icmd]+1],"overdamped") ) {
      sim_type = 2; // hi-friction limit for kinetics
    }
  } else if( !strcmp(opt[opt_ptr[icmd]],"temp") ) { // set the temperature
    set_temp(atof(opt[opt_ptr[icmd]+1]));

  } else if( !strcmp(opt[opt_ptr[icmd]],"nstep") ) { // # of steps
    nstep = atof(opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"istep_restart") ) { // where to restart from
    istep_restart = atof(opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"nup") ) { // # of steps before an update
    nup = atoi(opt[opt_ptr[icmd]+1]);

  } else if( !strcmp(opt[opt_ptr[icmd]],"run") ) { // set current run
    run = atoi((opt[opt_ptr[icmd]+1]));
    generator.set_seed(-100-run);

  } else if( !strcmp(opt[opt_ptr[icmd]],"ufname") ) { // set update file name
    strcpy(ufname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"rcfname") ) { // set restart coordinate file name
    strcpy(rcfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"cfname") ) { // set save coordinate file name
    strcpy(cfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"rgenfname") ) { // set random generator file name
    generator.set_fname(opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"unccfname") ) { // set save coordinate file name
    strcpy(unccfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"vfname") ) { // set save velocity file name
    strcpy(vfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"binfname") ) { // set save trajectory file name
    strcpy(binfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"uncbinfname") ) { // set save trajectory file name
    strcpy(uncbinfname,opt[opt_ptr[icmd]+1]);

  } else if( !strcmp(opt[opt_ptr[icmd]],"cutofftype") ) { // neighbor list on or off?
    if( !strcmp(opt[opt_ptr[icmd]+1],"neighborlist" ) ) { neighborlist = 1; }
    else { }

  } else if( !strcmp(opt[opt_ptr[icmd]],"nnlup") ) { // neighbor / cell list update frequency
    nnlup = atoi(opt[opt_ptr[icmd]+1]);

  } else if( !strcmp(opt[opt_ptr[icmd]],"boxl") ) { // box length for pbc                                            
    boxl = atof(opt[opt_ptr[icmd]+1]);

  } else if( !strcmp(opt[opt_ptr[icmd]],"restart") ) { // restart on or off?
    if( !strcmp(opt[opt_ptr[icmd]+1],"on" ) ) { restart = 1; }
    else { restart = 0; }
  } else if( !strcmp(opt[opt_ptr[icmd]],"rgen_restart") ) { // restart the generator?
    if( !strcmp(opt[opt_ptr[icmd]+1],"on" ) ) { rgen_restart = 1; }
    else { rgen_restart = 0; }
  } else if( !strcmp(opt[opt_ptr[icmd]],"t_step") ) {
    h = atof((opt[opt_ptr[icmd]+1]));
    halfh = h/2.0;
    a1 = h*(1.0-zeta*halfh);
    a2 = h*halfh;
    a3 = (1.0-h*zeta/2.0+(h*zeta)*(h*zeta)/4.0)/h;
    a4 = halfh*(1.0-h*zeta/2.0);
    a5 = h/zeta;
  } else if( !strcmp(opt[opt_ptr[icmd]],"zeta") ) { // friction coefficient
    zeta = atof((opt[opt_ptr[icmd]+1]));
    a1 = h*(1.0-zeta*halfh);
    a3 = (1.0-h*zeta/2.0+(h*zeta)*(h*zeta)/4.0)/h;
    a4 = halfh*(1.0-h*zeta/2.0);
    a5 = h/zeta;
  } else if( !strcmp(opt[opt_ptr[icmd]],"ion_conc") ) { // 
     ion_conc = atof((opt[opt_ptr[icmd]+1]));
     ion_str = ((ion_conc*z_p*z_p*Mol)+
		(ion_conc*z_n*z_n*Mol))*elec*elec/2.0;
     debye = sqrt(epsilon_r*T*KCpmtoJ/(8.0*pi*kelec*ion_str))*meter2angstrom;
     if( debye < Debye_limit ) debye = Debye_limit;
  } else {};
  
}

void set_temp(double temp)
{
  using namespace std;

  T = temp;
  debye = sqrt(epsilon_r*T*KCpmtoJ/(8.0*pi*kelec*ion_str))*meter2angstrom;
  if( stcks_allocated ) {
    for( int i=1; i<=nstck; i++ ) {
      rna_stcks[i].delta_G = rna_stcks[i].delta_H - 
	T * kcalpmol2K * rna_stcks[i].delta_S;
    }
  }
  
}

