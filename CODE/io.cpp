#include <iostream>
#include <fstream>
#include "io.h"
#include "global.h"
#include "param.h"

void read_input(const char* const ifile)
{

  using namespace std;

  ifstream in;
  char line[1024];
  char* tokPtr;
  char term = ';'; // terminates a command
  int newcmd = 1;
  int icmd;
  int nopt_tot = 0;
  int iopt;

  ncmd = 0;
  in.clear();
  in.open(ifile,ios::in);
  while(1) {
    in.getline(line,1024);
    if( in.eof() ) break;
    tokPtr = strtok(line," ");
    if( strchr(tokPtr,term)!=NULL ) { ncmd++; }
    while( tokPtr = strtok(NULL," ") ) {
      if( strchr(tokPtr,term)!=NULL ) { ncmd++; }
    }      
  }
  in.close();

  in.clear();
  in.open(ifile,ios::in);
  icmd = 0;
  while(1) {
    in.getline(line,1024);
    if( in.eof() ) break;
    tokPtr = strtok(line," ");;
    if( newcmd ) { 
      icmd++; 
      strcpy( cmd[icmd],tokPtr ); 
      opt_ptr[icmd] = nopt_tot+1; 
      newcmd = 0;
    } else { 
      nopt_tot++; 
      strcpy( opt[nopt_tot],tokPtr ); 
    }
    if( strchr(tokPtr,term)!=NULL ) {
      newcmd = 1;
    }
    while( tokPtr = strtok(NULL," ") ) {
      if( newcmd ) { 
        icmd++; 
        strcpy( cmd[icmd],tokPtr ); 
        opt_ptr[icmd] = nopt_tot+1; 
        newcmd = 0;
      } else { 
        nopt_tot++; 
        strcpy( opt[nopt_tot],tokPtr ); 
      }
      if( strchr(tokPtr,term)!=NULL ) {
        newcmd = 1;
      }
    }
  }
  opt_ptr[ncmd+1] = nopt_tot + 1;
  in.close();

  for( int icmd=1; icmd<=ncmd; icmd++ ) {
    for( int i=0; i<strlen(cmd[icmd]); i++ ) {
      if( cmd[icmd][i] == term ) cmd[icmd][i] = '\0';
    }

    for( int iopt = opt_ptr[icmd]; iopt < opt_ptr[icmd+1]; iopt++ ) {
      for( int i=0; i<strlen(opt[iopt]); i++ ) {
        if( opt[iopt][i] == ';' ) opt[iopt][i] = '\0';
      }

    }
  }
  
}

void load(int icmd)
{
  
  using namespace std;

  ifstream in;
  char line[2048];
  char* tokPtr;
  int test;
  int test1,test2;
  int ncon_tot;
  int icon_att, icon_rep;
  int i,j,k,l;
  int ibead,jbead;
  double real_phi, ideal_phi;
  double r_ij;
  int istate;

  if( !strcmp(opt[opt_ptr[icmd]],"bonds") ) { // load bonds
    cout << "[Reading in bonds...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in); // open file
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    nbnd = atoi(tokPtr); // read in number of bonds
    init_bonds(nbnd);
    for( int i=1; i<=nbnd; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      ibead_bnd[i] = atoi(tokPtr); // first bead index
      tokPtr = strtok(NULL," ");
      jbead_bnd[i] = atoi(tokPtr); // second bead index
      tokPtr = strtok(NULL," ");
      pdb_dist[i] = atof(tokPtr); // equilibrium distance (angstrom)
    }
    in.close(); // close file
    cout << "[Finished reading bonds (" << nbnd <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"angles")) { // load angles
    cout << "[Reading in angles...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    nang = atoi(tokPtr); // read in number of angles
    init_angles(nang);
    for( int i=1; i<=nang; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      ibead_ang[i] = atoi(tokPtr); // first bead index
      tokPtr = strtok(NULL," ");
      jbead_ang[i] = atoi(tokPtr); // second bead index
      tokPtr = strtok(NULL," ");
      kbead_ang[i] = atoi(tokPtr); // third bead index
      tokPtr = strtok(NULL," ");
      pdb_ang[i] = atof(tokPtr); // equilibrium angle (radians)
    }
    in.close();    
    cout << "[Finished reading angles (" << nang <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"torsions")) { // load torsions
    cout << "[Reading in torsions...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    ntor = atoi(tokPtr); // read in number of torsions
    release_torsions();
    init_torsions(ntor);
    for( int i=1; i<=ntor; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      ibead_tor[i] = atoi(tokPtr); // first bead index
      tokPtr = strtok(NULL," ");
      jbead_tor[i] = atoi(tokPtr); // second bead index
      tokPtr = strtok(NULL," ");
      kbead_tor[i] = atoi(tokPtr); // third bead index
      tokPtr = strtok(NULL," ");
      lbead_tor[i] = atoi(tokPtr); // fourth bead index
      tokPtr = strtok(NULL," ");
      real_phi = atof(tokPtr); // equilibrium torsion (angstrom)

      istate=(int)(3.0*real_phi/(2.0*pi))+1;
      if( istate==1 ) {
	ideal_phi = 60.5*pi/180.0;
      } else if( istate==2 ) {
	ideal_phi = 180.0*pi/180.0;
      } else if( istate==3 ) {
	ideal_phi = 299.5*pi/180.0;
      }
      tokPtr = strtok(NULL," ");
      cos_dih_ideal_dev[i] = cos(real_phi-ideal_phi);
      sin_dih_ideal_dev[i] = sin(real_phi-ideal_phi);

      tokPtr = strtok(NULL," ");
      alpha[i] = atof(tokPtr);
      tokPtr = strtok(NULL," ");
      beta[i] = atof(tokPtr);
      tokPtr = strtok(NULL," ");
      gam[i] = atof(tokPtr);
      tokPtr = strtok(NULL," ");
      delta[i] = atof(tokPtr);
    }
    in.close();
    cout << "[Finished reading torsions (" << ntor <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"vdw")) { // load rna-rna vdw
    cout << "[Reading in VDW interactions...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    ncon_att = atoi(tokPtr);
    tokPtr = strtok(NULL," ");
    tokPtr = strtok(NULL," ");
    ncon_rep = atoi(tokPtr);
    release_lj();
    init_lj(ncon_att, ncon_rep);
    ncon_tot = ncon_att + ncon_rep;
    icon_att = 0;
    icon_rep = 0;
    for( int i=1; i<=ncon_tot; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      ibead = atoi(tokPtr);
      tokPtr = strtok(NULL," ");
      jbead = atoi(tokPtr);
      tokPtr = strtok(NULL," ");
      r_ij = atof(tokPtr);
      if( (r_ij<rcut_nat)&&(is_base(ibead))&&(is_base(jbead)) ) {
	icon_att++;
	ibead_lj_nat[icon_att] = ibead;
	jbead_lj_nat[icon_att] = jbead;
	lj_nat_pdb_dist[icon_att] = r_ij;
	lj_nat_pdb_dist2[icon_att] = r_ij*r_ij;
	lj_nat_pdb_dist6[icon_att] = r_ij*r_ij*r_ij*r_ij*r_ij*r_ij;
	lj_nat_pdb_dist12[icon_att] = lj_nat_pdb_dist6[icon_att]*
	  lj_nat_pdb_dist6[icon_att];

	nil_att++;
	ibead_pair_list_att[nil_att] = ibead;
	jbead_pair_list_att[nil_att] = jbead;
	pl_lj_nat_pdb_dist[nil_att] = r_ij;
	pl_lj_nat_pdb_dist2[nil_att] = lj_nat_pdb_dist2[icon_att];
	pl_lj_nat_pdb_dist6[nil_att] = lj_nat_pdb_dist6[icon_att];
	pl_lj_nat_pdb_dist12[nil_att] = lj_nat_pdb_dist12[icon_att];
      } else if(is_base(ibead)&&is_base(jbead)) {
	icon_rep++;
	ibead_lj_non_nat[icon_rep] = ibead;
	jbead_lj_non_nat[icon_rep] = jbead;

	nil_rep++;
	ibead_pair_list_rep[nil_rep] = ibead;
	jbead_pair_list_rep[nil_rep] = jbead;

	switch_fnb[icon_rep] = 1;
      } else {
	icon_rep++;
	ibead_lj_non_nat[icon_rep] = ibead;
	jbead_lj_non_nat[icon_rep] = jbead;

	nil_rep++;
	ibead_pair_list_rep[nil_rep] = ibead;
	jbead_pair_list_rep[nil_rep] = jbead;

	switch_fnb[icon_rep] = 0;
      }
    }
    in.close();
    cout << "[Finished reading VDW interactions (" << icon_att << "/" << icon_rep <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"stacks")) { // load stacks
    cout << "[Reading in stacking...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    nstck = atoi(tokPtr);
    release_stacks();
    init_stacks(nstck);
    for( int i=1; i<=nstck; i++ ) {
      for( int j=1; j<=4; j++ ) {
	in.getline(line,2048);
	tokPtr = strtok(line," ");
	rna_stcks[i].ibead_ang[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].jbead_ang[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].kbead_ang[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].pdb_ang[j] = atof(tokPtr);
	rna_stcks[i].cos_pdb_ang[j] = cos(rna_stcks[i].pdb_ang[j]);
	rna_stcks[i].sin_pdb_ang[j] = sin(rna_stcks[i].pdb_ang[j]);
      }
      for( int j=1; j<=2; j++ ) {
	in.getline(line,2048);
	tokPtr = strtok(line," ");
	rna_stcks[i].ibead_dist[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].jbead_dist[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].pdb_dist[j] = atof(tokPtr);
      }      
      for( int j=1; j<=2; j++ ) {
	in.getline(line,2048);
	tokPtr = strtok(line," ");
	rna_stcks[i].ibead_tor[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].jbead_tor[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].kbead_tor[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].lbead_tor[j] = atoi(tokPtr);
	tokPtr = strtok(NULL," ");
	rna_stcks[i].pdb_tor[j] = atof(tokPtr);
	rna_stcks[i].cos_pdb_tor[j] = cos(rna_stcks[i].pdb_tor[j]);
	rna_stcks[i].sin_pdb_tor[j] = sin(rna_stcks[i].pdb_tor[j]);
      }
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      rna_stcks[i].delta_H = atof(tokPtr);
      tokPtr = strtok(NULL," ");
      rna_stcks[i].delta_S = atof(tokPtr);
      rna_stcks[i].delta_G = rna_stcks[i].delta_H - 
	T * kcalpmol2K * rna_stcks[i].delta_S;
    }
    in.close();
    cout << "[Finished reading stacking (" << nstck <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"elec")) { // load electrostatics
    cout << "[Reading in electrostatics...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    nelec = atoi(tokPtr);
    release_elec();
    init_elec(nelec);
    for( int i=1; i<=nelec; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      ibead_elec[i] = atoi(tokPtr);
      tokPtr = strtok(NULL," ");
      jbead_elec[i] = atoi(tokPtr);
    }
    in.close();
    cout << "[Finished reading electrostatics (" << nelec <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"init")) { // load init coordinates
    cout << "[Reading in initial coordinates...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    nbead = atoi(tokPtr); // read in number of beads
    init_pos(nbead);
    for( int i=1; i<=nbead; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      tokPtr = strtok(NULL," ");
      pos[i].x = atof(tokPtr);
      unc_pos[i].x = pos[i].x;
      tokPtr = strtok(NULL," ");
      pos[i].y = atof(tokPtr);
      unc_pos[i].y = pos[i].y;
      tokPtr = strtok(NULL," ");
      pos[i].z = atof(tokPtr);
      unc_pos[i].z = pos[i].z;
    }
    in.close();
    cout << "[Finished reading initial coordinates (" << nbead << ")]" << endl;
  }
  
}

void record_traj(char* fname,char* fname2)
{

  using namespace std;

  char oline[1024];
  char oline2[1024];
  ofstream trajfile;
  ofstream trajfile2;

  trajfile.open(fname,ios::out | ios::app);
  trajfile2.open(fname2,ios::out | ios::app);

  for( int i=1; i<=nbead; i++ ) {
    sprintf(oline,"%f %f %f",pos[i].x,pos[i].y,pos[i].z);
    sprintf(oline2,"%f %f %f",unc_pos[i].x,unc_pos[i].y,unc_pos[i].z);
    trajfile << oline << endl;
    trajfile2 << oline2 << endl;

  }

  trajfile.close();
  trajfile2.close();

} 

void save_coords(char* fname,char* fname2)
{

  using namespace std;

  char oline[1024];
  ofstream ofile;
  ofstream ofile2;

  ofile.open(fname,ios::out);
  ofile2.open(fname2,ios::out);
  for( int i=1; i<=nbead; i++ ) {
    sprintf(oline,"%d %f %f %f",i,pos[i].x,
	    pos[i].y,pos[i].z);
    ofile << oline << endl;
    sprintf(oline,"%d %f %f %f",i,unc_pos[i].x,
	    unc_pos[i].y,unc_pos[i].z);
    ofile2 << oline << endl;
  }
  ofile.close();
  ofile2.close();
  
}

void save_unccoords(char* fname)
{

  using namespace std;

  char oline[1024];
  ofstream ofile;

  ofile.open(fname,ios::out);
  for( int i=1; i<=nbead; i++ ) {
    sprintf(oline,"%d %f %f %f",i,unc_pos[i].x,
	    unc_pos[i].y,unc_pos[i].z);
    ofile << oline << endl;
  }
  ofile.close();

}

void load_coords(char* fname,char* fname2)
{

  using namespace std;

  char iline[1024];
  ifstream ifile;
  ifstream ifile2;
  char* tokPtr;

  ifile.clear();
  ifile2.clear();
  ifile.open(fname,ios::in);
  ifile2.open(fname2,ios::in);
  for( int i=1; i<=nbead; i++ ) {
    ifile.getline(iline,1024);
    tokPtr = strtok(iline," ");
    tokPtr = strtok(NULL," ");
    pos[i].x = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    pos[i].y = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    pos[i].z = atof(tokPtr);
    ifile2.getline(iline,1024);
    tokPtr = strtok(iline," ");
    tokPtr = strtok(NULL," ");
    unc_pos[i].x = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    unc_pos[i].y = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    unc_pos[i].z = atof(tokPtr);
  }
  ifile.close();
  ifile2.close();

}

void save_vels(char* fname)
{

  using namespace std;

  char oline[1024];
  ofstream ofile;

  ofile.open(fname,ios::out);
  for( int i=1; i<=nbead; i++ ) {
    sprintf(oline,"%d %f %f %f",i,vel[i].x,
	    vel[i].y,vel[i].z);
    ofile << oline << endl;
  }
  ofile.close();

}

void load_vels(char* fname)
{

  using namespace std;

  char iline[1024];
  ifstream ifile;
  char* tokPtr;

  ifile.clear();
  ifile.open(fname,ios::in);
  for( int i=1; i<=nbead; i++ ) {
    ifile.getline(iline,1024);
    tokPtr = strtok(iline," ");
    tokPtr = strtok(NULL," ");
    vel[i].x = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    vel[i].y = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    vel[i].z = atof(tokPtr);

  }
  ifile.close();

}

int is_base(int baseno)
{

  using namespace std;

  int retval = 0;

  for( int i=2; i<=nbead; i+=3 ) {
    if( i==baseno ) {
      retval=1; break;
    }
  }

  return retval;

}

int is_phosphate(int phosno)
{
  
  using namespace std;

  int retval = 0;
  
  for( int i=3; i<=nbead; i+=3 ) {
    if( i==phosno ) {
      retval=1; break;
    }
  }
  
  return retval;
  
}

