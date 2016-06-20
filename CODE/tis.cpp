 #include <iostream>
 #include <fstream>
 #include <cstring>
 #include <math.h>
 #include <cstdlib>
 #include <ctime>
 #include <cstdio>
 #include <unistd.h>
 #include "tis.h"
 #include "energy.h"
 #include "global.h"
 #include "io.h"
 #include "param.h"
 #include "misc.h"
 #include "random_generator.h"

 int main(int argc,char* argv[])
 {

   using namespace std;

   if( argc<2 ) {
     cerr << "Usage: " << argv[0] <<  " < input_file >" << endl;
     exit(-1);    
   }
   time_t tm0 = time(0); // wall time at this point
   clock_t ck0 = clock(); // clock ticks at this point
   cout << "CURRENT TIME IS: " << ctime(&tm0);
   if( getcwd(pathname,MAXPATHLEN)==NULL ) {
     cerr << "PROBLEM GETTING PATH" << endl;
   } else {
     cout << "CURRENT WORKING DIRECTORY: " << pathname << endl;
   }

   alloc_arrays(); // allocates certain arrays and initializes some variables
   read_input(argv[1]); // read input file

   ex_cmds(); // perform commands (simulation)

   time_t tm1 = time(0);
   clock_t ck1 = clock();
   cout << "+-------------------+" << endl;
   cout << "| Simulation Stats: |" << endl;
   cout << "+-------------------+" << endl;
   cout << "Wall Time              : " << difftime(tm1,tm0) << " sec" << endl;
   cout << "Total Computation Time : " << float(ck1-ck0)/CLOCKS_PER_SEC << " sec" << endl;
   cout << "Computation Rate       : " << float(ck1-ck0)/CLOCKS_PER_SEC/nstep << " sec / timestep" << endl;
   cout << "CURRENT TIME IS        : " << ctime(&tm1);

   return 0;

 }

 void ex_cmds()
 {

   using namespace std;

   char oline[1024];
   int iopt;

   for( int i=1; i<=ncmd; i++ ) {
      // read data
      if( !strcmp(cmd[i],"load") ) { load(i); }
      // set parameters
      else if( !strcmp(cmd[i],"set") ) { set_params(i); }
      // run simulation
      else if( !strcmp(cmd[i],"run") ) { simulation_ctrl(); }
      // ???
      else {};
   }

 }

 void simulation_ctrl()
 {

   using namespace std;

   switch( sim_type ) {
   case 1:
     underdamped_ctrl();
     break;
   case 2:
     overdamped_ctrl();
     break;
   default:
     cerr << "UNRECOGNIZED SIM_TYPE!" << endl;
     exit(-1);
   }
 }

 void underdamped_ctrl()
 {

   using namespace std;

   char oline[2048];
   double istep = 1.0;
   int iup = 1;
   ofstream out(ufname,ios::out|ios::app);
   static int first_time = 1;

   coord* incr = new coord[nbead+1];

   if( (!restart)&&first_time ) { // zero out the velocities and forces
     for( int i=1; i<=nbead; i++ ) {
       vel[i].x = 0.0;
       vel[i].y = 0.0;
       vel[i].z = 0.0;
       force[i].x = 0.0;
       force[i].y = 0.0;
       force[i].z = 0.0;
     }
   }

   print_sim_params();

   if (neighborlist == 1) {
     update_neighbor_list();
     update_pair_list();
   }

   set_potential();
   set_forces();

   char line[2048];

   if( restart ) {
     load_coords(cfname,unccfname);
     load_vels(vfname);
     istep = istep_restart + 1.0;
   }

   if( rgen_restart ) {
     generator.restart();
   }

   if( first_time ) {

     energy_eval();
     force_eval();

   }

   if( binsave ) {
     if( (first_time)&&(!rgen_restart) ) {
       record_traj(binfname,uncbinfname);     
     }
  
     while( istep <= nstep ) {

       if ((inlup % nnlup) == 0) {
	 if (neighborlist == 1) {
	   update_neighbor_list();
	 }
	 //	 fprintf(stderr, "(%.0lf) neighbor list: (%d/%d)\n", istep, nnl_att, nnl_rep);
	 inlup = 0;
       }
       inlup++;

       if (neighborlist == 1) {
	 update_pair_list();
	 //	 fprintf(stderr, "(%.0lf) pair list: (%d/%d)\n", istep, nil_att, nil_rep);
       }

       underdamped_iteration(incr);
       if( !(iup%nup) ) { // updates

	 energy_eval();
	 calculate_observables(incr);
	 sprintf(oline,"%.0lf %f %f %f %f %f %f %f %f %f %f %d %f",
		 istep,T,kinT,e_bnd,e_ang,e_tor,e_stack,e_vdw_rr,e_elec,rna_etot,
		 Q,contct_nat,rgsq);
	 out << oline << endl;
	 iup = 0;
	 record_traj(binfname,uncbinfname);     
	 save_coords(cfname,unccfname);
	 save_vels(vfname);
	 generator.save_state();

       }
       istep += 1.0;
       iup++;
      
    }
    out.close();
  }
  
  if( first_time ) first_time = 0;
  
  delete [] incr;

  return;

}

void overdamped_ctrl()
{

  using namespace std;

  char oline[2048];
  double istep = 1.0;
  int iup = 1;
  ofstream out(ufname,ios::out|ios::app);
  static int first_time = 1;

  coord* incr = new coord[nbead+1];

  if( (!restart)&&first_time ) { // zero out the velocities and forces
    for( int i=1; i<=nbead; i++ ) {
      vel[i].x = 0.0;
      vel[i].y = 0.0;
      vel[i].z = 0.0;
      force[i].x = 0.0;
      force[i].y = 0.0;
      force[i].z = 0.0;
    }
  }
  set_potential();
  set_forces();
  
  sprintf(oline,"CURRENT TEMPERATURE: %f",T);
  cout << oline << endl;

  char line[2048];

  if( restart ) {
    load_coords(cfname,unccfname);
    load_vels(vfname);
    istep = istep_restart + 1.0;
  }
  
  if( rgen_restart ) {
    generator.restart();
  }

  if( first_time ) {

    energy_eval();
    force_eval();

  }

  if( binsave ) {
    if( (first_time)&&(!rgen_restart) ) {
      record_traj(binfname,uncbinfname);     
    }
    while( istep <= nstep ) {
      
      overdamped_iteration(incr);
      if( !(iup%nup) ) { // updates
	energy_eval();
	calculate_observables(incr);
        sprintf(oline,"%.0lf %f %f %f %f %f %f %f %f %f %f %d %f",
                istep,T,kinT,e_bnd,e_ang,e_tor,e_stack,e_vdw_rr,e_elec,rna_etot,
                Q,contct_nat,rgsq);
	out << oline << endl;
	iup = 0;
	record_traj(binfname,uncbinfname);     
	save_coords(cfname,unccfname);
	save_vels(vfname);
	generator.save_state();
      }
      istep += 1.0;
      iup++;
      
    }
    out.close();
  }
  
  if( first_time ) first_time = 0;
  
  delete [] incr;

  return;

}

void underdamped_iteration(coord* increment)
{
  
  using namespace std;
  
  static const double eps = 1.0e-5;

  for( int i=1; i<=nbead; i++ ) {
    
    // compute position increments
    
    increment[i].x = a1*vel[i].x + a2*force[i].x;
    increment[i].y = a1*vel[i].y + a2*force[i].y;
    increment[i].z = a1*vel[i].z + a2*force[i].z;

    // update bead positions

    pos[i].x += increment[i].x;
    pos[i].y += increment[i].y;
    pos[i].z += increment[i].z;

    pos[i].x -= boxl*rnd(pos[i].x/boxl);
    pos[i].y -= boxl*rnd(pos[i].y/boxl);
    pos[i].z -= boxl*rnd(pos[i].z/boxl);

    unc_pos[i].x += increment[i].x;
    unc_pos[i].y += increment[i].y;
    unc_pos[i].z += increment[i].z;

  }

  
  // force_update

  force_eval();

  if( T < eps ) return; // don't update velocities for steepest descent

  // update_velocities

  for( int i=1; i<=nbead; i++ ) {

    // compute velocity increments
    
    vel[i].x = a3*increment[i].x + a4*force[i].x;
    vel[i].y = a3*increment[i].y + a4*force[i].y;
    vel[i].z = a3*increment[i].z + a4*force[i].z;
    
  }
  
}

void overdamped_iteration(coord* increment)
{

  using namespace std;
  
  for( int i=1; i<=nbead; i++ ) {
    
    // compute position increments
    
    increment[i].x = a5*force[i].x;
    increment[i].y = a5*force[i].y;
    increment[i].z = a5*force[i].z;

    // update bead positions

    unc_pos[i].x += increment[i].x;
    unc_pos[i].y += increment[i].y;
    unc_pos[i].z += increment[i].z;

    pos[i].x += increment[i].x;
    pos[i].y += increment[i].y;
    pos[i].z += increment[i].z;

    pos[i].x -= boxl*rnd(pos[i].x/boxl);
    pos[i].y -= boxl*rnd(pos[i].y/boxl);
    pos[i].z -= boxl*rnd(pos[i].z/boxl);

  }

  // force_update

  force_eval();

}



void calculate_observables(coord* increment)
{

  using namespace std;

  char oline[1024];
  double dx,dy,dz,d;
  static const double tol = 1.0; // tolerance for chi distances
  static const double chinorm = (double(nbead*nbead)-5.0*double(nbead)+6.0)/2.0;
  double sumvsq;
  int nchi;
  int ibead, jbead;
  float r_ij;
  char line[2048];

  // chi, contct_nat, contct_tot, Q
  
  contct_nat = 0;
  for( int i=1; i<=ncon_att; i++ ) {

    ibead = ibead_lj_nat[i];
    jbead = jbead_lj_nat[i];
    r_ij = lj_nat_pdb_dist[i];

    dx = unc_pos[ibead].x-unc_pos[jbead].x;
    dy = unc_pos[ibead].y-unc_pos[jbead].y;
    dz = unc_pos[ibead].z-unc_pos[jbead].z;

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d = sqrt( dx*dx+dy*dy+dz*dz );
    if(d/r_ij < 1.25) {
      contct_nat++;
    }
  }
  Q = double(contct_nat)/ncon_att;

  // rgsq
  
  rgsq = 0.0;
  for( int i=1; i<=nbead-1; i++ ) {
    for( int j=i+1; j<=nbead; j++ ) {
      dx = unc_pos[i].x-unc_pos[j].x;
      dy = unc_pos[i].y-unc_pos[j].y;
      dz = unc_pos[i].z-unc_pos[j].z;

      dx -= boxl*rnd(dx/boxl);
      dy -= boxl*rnd(dy/boxl);
      dz -= boxl*rnd(dz/boxl);

      rgsq += (dx*dx+dy*dy+dz*dz);
    }
  }
  rgsq /= double(nbead*nbead);
  
  // kinT

  if( sim_type == 1 ) {
    sumvsq = 0.0;
    for( int i=1; i<=nbead; i++ ) {
      sumvsq += vel[i].x*vel[i].x
	+ vel[i].y*vel[i].y
	+ vel[i].z*vel[i].z;
    }
    kinT = sumvsq/(3.0*double(nbead));
  } else if( sim_type == 2 ) {
    sumvsq = 0.0;
    for( int i=1; i<=nbead; i++ ) {
      sumvsq += increment[i].x*increment[i].x +
	increment[i].y*increment[i].y +
	increment[i].z*increment[i].z;
    }
    sumvsq *= zeta/(2.0*h);
    kinT = sumvsq/(3.0*double(nbead));
  } else {}
  
  
}

void print_sim_params() {

  using namespace std;

  char oline[2048];

  cout << endl;
  sprintf(oline,"+------------------------+");
  cout << oline << endl;
  sprintf(oline,"| Simulation Parameters: |");
  cout << oline << endl;
  sprintf(oline,"+------------------------+");
  cout << oline << endl;

  if (sim_type == 1) {
    sprintf(oline,"Simulation Type                   : %s", "Underdamped");
    cout << oline << endl;
  } else if (sim_type == 2) {
    sprintf(oline,"Simulation Type                   : %s", "Overdamped");
    cout << oline << endl;
  } else {
    cerr << "UNRECOGNIZED SIMULATION TYPE!" << endl;
    exit(-1);
  }

  sprintf(oline,"Simulation Temperature            : %.3f",T);
  cout << oline << endl;

  sprintf(oline,"Start Time Step                   : %.0lf", istep_restart);
  cout << oline << endl;

  sprintf(oline,"Final Time Step                   : %.0lf", nstep);
  cout << oline << endl;

  sprintf(oline,"Output Frequency                  : %d", nup);
  cout << oline << endl;

  sprintf(oline,"Friction Coefficient              : %.0e", zeta);
  cout << oline << endl;

  sprintf(oline,"PBC Box Length                    : %.1f", boxl);
  cout << oline << endl;

  if (neighborlist == 1) {
    sprintf(oline,"Long-range Cutoff Type            : %s", "Neighbor List");
    cout << oline << endl;
    sprintf(oline,"Neighbor List Update Frequency    : %d", nnlup);
    cout << oline << endl;
  } else {
    sprintf(oline,"Long-range Cutoff Type            : %s", "None");
    cout << oline << endl;
  }

  cout << endl;

}

void update_neighbor_list() {

  double dx, dy, dz;
  double d2;
  int ibead, jbead;
  double rcut, rcut2;

  nnl_att = 0;
  nnl_rep = 0;

  for (int i=1; i<=ncon_att; i++) {

    ibead = ibead_lj_nat[i];
    jbead = jbead_lj_nat[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    d2 = dx*dx+dy*dy+dz*dz;

    rcut = 3.2*lj_nat_pdb_dist[i];
    rcut2 = rcut*rcut;

    if (d2 < rcut2) {
      // add to neighbor list
      nnl_att++;
      ibead_neighbor_list_att[nnl_att] = ibead;
      jbead_neighbor_list_att[nnl_att] = jbead;
      nl_lj_nat_pdb_dist[nnl_att] = lj_nat_pdb_dist[i];
      nl_lj_nat_pdb_dist2[nnl_att] = lj_nat_pdb_dist2[i];
      nl_lj_nat_pdb_dist6[nnl_att] = lj_nat_pdb_dist6[i];
      nl_lj_nat_pdb_dist12[nnl_att] = lj_nat_pdb_dist12[i];
    }
  }

  for (int i=1; i<=ncon_rep; i++) {

    ibead = ibead_lj_non_nat[i];
    jbead = jbead_lj_non_nat[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    d2 = dx*dx+dy*dy+dz*dz;

    rcut = 3.2*sigma_rep;
    rcut2 = rcut*rcut;

    if (d2 < rcut2) {
      // add to neighbor list
      nnl_rep++;
      ibead_neighbor_list_rep[nnl_rep] = ibead;
      jbead_neighbor_list_rep[nnl_rep] = jbead;
    }

  }
}

void update_pair_list() {

  using namespace std;

  // declare host variables
  double dx, dy, dz;
  double d2;
  unsigned int ibead, jbead;
  double rcut, rcut2;

  nil_att = 0;
  nil_rep = 0;

  // declare device variables

  // should be native distance
  for (int i=1; i<=nnl_att; i++) {

    ibead = ibead_neighbor_list_att[i];
    jbead = jbead_neighbor_list_att[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    d2 = dx*dx+dy*dy+dz*dz;

    rcut = 2.5*nl_lj_nat_pdb_dist[i];
    rcut2 = rcut*rcut;

    if (d2 < rcut2) {
      // add to interaction pair list
      nil_att++;
      ibead_pair_list_att[nil_att] = ibead;
      jbead_pair_list_att[nil_att] = jbead;
      pl_lj_nat_pdb_dist[nil_att] = nl_lj_nat_pdb_dist[i];
      pl_lj_nat_pdb_dist2[nil_att] = nl_lj_nat_pdb_dist2[i];
      pl_lj_nat_pdb_dist6[nil_att] = nl_lj_nat_pdb_dist6[i];
      pl_lj_nat_pdb_dist12[nil_att] = nl_lj_nat_pdb_dist12[i];
    }

  }

  for (int i=1; i<=nnl_rep; i++) {

    ibead = ibead_neighbor_list_rep[i];
    jbead = jbead_neighbor_list_rep[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    d2 = dx*dx+dy*dy+dz*dz;

    rcut = 2.5*sigma_rep;
    rcut2 = rcut*rcut;

    if (d2 < rcut2) {
      // add to interaction pair list
      nil_rep++;
      ibead_pair_list_rep[nil_rep] = ibead;
      jbead_pair_list_rep[nil_rep] = jbead;
    }
  }
}

