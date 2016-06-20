#include <cmath>
#include <cstring>
#include "energy.h"
#include "global.h"
#include "tis.h"
#include "random_generator.h"
#include "misc.h"

void set_potential() {

  using namespace std;
  
  int iterm;
  
  iterm = 0;
  for( int i=1; i<=mpot_term; i++ ) {
    switch(i) {
    case 1:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &bond_energy;
      }
      break;
    case 2:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &angular_energy;
      }
      break;
    case 3:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &torsional_energy;
      }
      break;
    case 4:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &stacking_energy;
      }
      break;
    case 5:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &vdw_energy;
      }
      break;
    case 6:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &electrostatic_energy;
      }
      break;
    default:
      break;
    }
  }

}

void set_forces()
{
  
  using namespace std;

  int iterm;

  iterm = 0;
  for( int i=1; i<=mforce_term; i++ ) {
    switch(i) {
    case 1:
      if( force_term_on[i] ) {
	force_term[++iterm] = &random_force;
      }
      break;
    case 2:
      if( force_term_on[i] ) {
	force_term[++iterm] = &bond_forces;
      } 
      break;
    case 3:
      if( force_term_on[i] ) {
	force_term[++iterm] = &angular_forces;
      }
      break;
    case 4:
      if( force_term_on[i] ) {
	force_term[++iterm] = &torsional_forces;
      }
      break;
    case 5:
      if( force_term_on[i] ) {
	force_term[++iterm] = &stacking_forces;
      }
      break;
    case 6:
      if( force_term_on[i] ) {
	force_term[++iterm] = &vdw_forces;
      }
      break;
    case 7:
      if( force_term_on[i] ) {
	force_term[++iterm] = &electrostatic_forces;
      }
      break;
    default:
      break;
    }
  }
  
}

void clear_forces() {
  
  using namespace std;

  for( int i=1; i<=nbead; i++ ) {
    force[i].x = 0.0;
    force[i].y = 0.0;
    force[i].z = 0.0;
  }
  
}

void energy_eval()
{
  
  using namespace std;
  char oline[1024];

  for( int i=1; i<=npot_term; i++ ) {
    pot_term[i]();
  }

  rna_etot = e_bnd + e_ang + e_tor + e_stack + e_vdw_rr + e_elec;
  system_etot = rna_etot;

}

void force_eval()
{

  using namespace std;
  char oline[1024];

  clear_forces();
  
  for( int i=1; i<=nforce_term; i++ ) {
    force_term[i]();
  }

}

void random_force() {

  using namespace std;
  
  double var;

  var = sqrt(2.0*T*zeta/h);
  
  for( int i=1; i<=nbead; i++ ) {
    force[i].x += var*generator.gasdev();
    force[i].y += var*generator.gasdev();
    force[i].z += var*generator.gasdev();

  }
  
}

void bond_energy()
{
  
  using namespace std;

  int ibead, jbead;
  double dx, dy, dz, d;
  char line[2048];

  e_bnd = 0.0;
  for( int i=1; i<=nbnd; i++ ) {
    
    ibead = ibead_bnd[i];
    jbead = jbead_bnd[i];

    dx = unc_pos[jbead].x-unc_pos[ibead].x;
    dy = unc_pos[jbead].y-unc_pos[ibead].y;
    dz = unc_pos[jbead].z-unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d = sqrt(dx*dx+dy*dy+dz*dz);
    
    e_bnd += (d-pdb_dist[i])*(d-pdb_dist[i]);

  }

  e_bnd *= e_bnd_coeff;

  return;

}

void bond_forces()
{

  using namespace std;

  int ibead, jbead;
  double dx, dy, dz, d;
  double fx, fy, fz;
  const static double tol = 1e-7;
  double co;

  char line[2048];

  for( int i=1; i<=nbnd; i++ ) {
    
    ibead = ibead_bnd[i];
    jbead = jbead_bnd[i];
    
    dx = unc_pos[jbead].x-unc_pos[ibead].x;
    dy = unc_pos[jbead].y-unc_pos[ibead].y;
    dz = unc_pos[jbead].z-unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d = sqrt(dx*dx+dy*dy+dz*dz);

    if( d<tol ) continue;
    co = k_bnd*(d-pdb_dist[i])/d;
    fx = co*dx;
    fy = co*dy;
    fz = co*dz;    

    force[ibead].x += fx;
    force[ibead].y += fy;
    force[ibead].z += fz;

    force[jbead].x -= fx;
    force[jbead].y -= fy;
    force[jbead].z -= fz;

  }
  
}

void angular_energy()
{

  using namespace std;

  e_ang = 0.0;
  int ibead, jbead, kbead;
  coord r_ij, r_jk;
  double C_12; // dot product
  double d1,d2; // distances
  double cos_theta,theta;

  for( int i=1; i<=nang; i++ ) {
    
    ibead = ibead_ang[i];
    jbead = jbead_ang[i];
    kbead = kbead_ang[i];
    
    r_ij.x = unc_pos[ibead].x - unc_pos[jbead].x;
    r_ij.y = unc_pos[ibead].y - unc_pos[jbead].y;
    r_ij.z = unc_pos[ibead].z - unc_pos[jbead].z;

    // min images

    r_ij.x -= boxl*rnd(r_ij.x/boxl);
    r_ij.y -= boxl*rnd(r_ij.y/boxl);
    r_ij.z -= boxl*rnd(r_ij.z/boxl);

    r_jk.x = unc_pos[kbead].x - unc_pos[jbead].x;
    r_jk.y = unc_pos[kbead].y - unc_pos[jbead].y;
    r_jk.z = unc_pos[kbead].z - unc_pos[jbead].z;

    // min images

    r_jk.x -= boxl*rnd(r_jk.x/boxl);
    r_jk.y -= boxl*rnd(r_jk.y/boxl);
    r_jk.z -= boxl*rnd(r_jk.z/boxl);

    d1 = sqrt(r_ij.x*r_ij.x + r_ij.y*r_ij.y + r_ij.z*r_ij.z);
    C_12 = r_ij.x*r_jk.x + r_ij.y*r_jk.y + r_ij.z*r_jk.z;
    d2 = sqrt(r_jk.x*r_jk.x + r_jk.y*r_jk.y + r_jk.z*r_jk.z);

    cos_theta = C_12/(d1*d2);
    theta = acos(cos_theta);

    e_ang += (theta-pdb_ang[i])*(theta-pdb_ang[i]);

  }

  e_ang *= e_ang_coeff;
  
  return;
  
}

void angular_forces()
{

  using namespace std;

  int ibead, jbead, kbead;
  coord r_ij, r_jk;
  double C_11,C_12,C_22; // dot products
  double d1,d2; // distances
  double cos_theta,sin_theta2,sin_theta,theta;
  double co1,co2,co3,co4,co5;
  const static double tol = 1e-7;

  char line[2048];
  
  for( int i=1; i<=nang; i++ ) {

    ibead = ibead_ang[i];
    jbead = jbead_ang[i];
    kbead = kbead_ang[i];
    
    r_ij.x = unc_pos[ibead].x - unc_pos[jbead].x;
    r_ij.y = unc_pos[ibead].y - unc_pos[jbead].y;
    r_ij.z = unc_pos[ibead].z - unc_pos[jbead].z;

    // min images

    r_ij.x -= boxl*rnd(r_ij.x/boxl);
    r_ij.y -= boxl*rnd(r_ij.y/boxl);
    r_ij.z -= boxl*rnd(r_ij.z/boxl);

    r_jk.x = unc_pos[kbead].x - unc_pos[jbead].x;
    r_jk.y = unc_pos[kbead].y - unc_pos[jbead].y;
    r_jk.z = unc_pos[kbead].z - unc_pos[jbead].z;

    // min images

    r_jk.x -= boxl*rnd(r_jk.x/boxl);
    r_jk.y -= boxl*rnd(r_jk.y/boxl);
    r_jk.z -= boxl*rnd(r_jk.z/boxl);

    C_11 = r_ij.x*r_ij.x + r_ij.y*r_ij.y + r_ij.z*r_ij.z;
    C_12 = r_ij.x*r_jk.x + r_ij.y*r_jk.y + r_ij.z*r_jk.z;
    C_22 = r_jk.x*r_jk.x + r_jk.y*r_jk.y + r_jk.z*r_jk.z;
    d1 = sqrt(C_11);
    d2 = sqrt(C_22);

    cos_theta = C_12/(d1*d2);
    sin_theta2 = 1 - cos_theta*cos_theta;
    if( sin_theta2 < tol ) continue;
    if( d1*d2 < tol ) continue;
    sin_theta = sqrt(sin_theta2);
    theta = acos(cos_theta);

    co1 = k_ang*(theta - pdb_ang[i])/sin_theta;
    co2 = 1.0/(d1*d2);
    co3 = co2*co2*co2;

    co4 = co3*C_12*C_22;
    force[ibead].x += co1*(r_jk.x*co2 - r_ij.x*co4);
    force[ibead].y += co1*(r_jk.y*co2 - r_ij.y*co4);
    force[ibead].z += co1*(r_jk.z*co2 - r_ij.z*co4);
    
    co4 = co3*C_12*C_11;
    force[kbead].x += co1*(r_ij.x*co2 - r_jk.x*co4);
    force[kbead].y += co1*(r_ij.y*co2 - r_jk.y*co4);
    force[kbead].z += co1*(r_ij.z*co2 - r_jk.z*co4);

    co4 = co3*C_12*C_22 - co2;
    co5 = co3*C_12*C_11 - co2;
    force[jbead].x += co1*(r_ij.x*co4 + r_jk.x*co5);
    force[jbead].y += co1*(r_ij.y*co4 + r_jk.y*co5);
    force[jbead].z += co1*(r_ij.z*co4 + r_jk.z*co5);

  }

}

void torsional_energy()
{

  using namespace std;

  int ibead,jbead,kbead,lbead;
  double dx1, dy1, dz1;
  double dx2, dy2, dz2;
  double dx3, dy3, dz3;
  double C_11,C_12,C_13;
  double C_22,C_23;
  double C_33;
  double D_12,D_23;
  double cosphi,sinphi,sign;
  double B;

  char line[2048];

  e_tor = 0.0;
  for( int itor=1; itor<=ntor; itor++ ) {
    
    ibead = ibead_tor[itor];
    jbead = jbead_tor[itor];
    kbead = kbead_tor[itor];
    lbead = lbead_tor[itor];
    
    dx1 = unc_pos[jbead].x - unc_pos[ibead].x;
    dy1 = unc_pos[jbead].y - unc_pos[ibead].y;
    dz1 = unc_pos[jbead].z - unc_pos[ibead].z;

    dx2 = unc_pos[kbead].x - unc_pos[jbead].x;
    dy2 = unc_pos[kbead].y - unc_pos[jbead].y;
    dz2 = unc_pos[kbead].z - unc_pos[jbead].z;
    
    dx3 = unc_pos[lbead].x - unc_pos[kbead].x;
    dy3 = unc_pos[lbead].y - unc_pos[kbead].y;
    dz3 = unc_pos[lbead].z - unc_pos[kbead].z;

    // min images

    dx1 -= boxl*rnd(dx1/boxl);
    dy1 -= boxl*rnd(dy1/boxl);
    dz1 -= boxl*rnd(dz1/boxl);
       
    dx2 -= boxl*rnd(dx2/boxl);
    dy2 -= boxl*rnd(dy2/boxl);
    dz2 -= boxl*rnd(dz2/boxl);
       
    dx3 -= boxl*rnd(dx3/boxl);
    dy3 -= boxl*rnd(dy3/boxl);
    dz3 -= boxl*rnd(dz3/boxl);
    
    C_11 = dx1*dx1 + dy1*dy1 + dz1*dz1;
    C_12 = dx1*dx2 + dy1*dy2 + dz1*dz2;
    C_13 = dx1*dx3 + dy1*dy3 + dz1*dz3;
    C_22 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    C_23 = dx2*dx3 + dy2*dy3 + dz2*dz3;
    C_33 = dx3*dx3 + dy3*dy3 + dz3*dz3;
    
    D_12 = C_11*C_22 - C_12*C_12;
    D_23 = C_22*C_33 - C_23*C_23;    

    cosphi = (C_12*C_23-C_13*C_22)/sqrt(D_12*D_23);
    sinphi = sqrt(1.0 - cosphi*cosphi);
    sign = (dz1*dy3-dy1*dz3)*dx2 + (dx1*dz3-dz1*dx3)*dy2
      + (dy1*dx3-dx1*dy3)*dz2;
    if( sign<0.0 ) sinphi = -sinphi;    

    B = cosphi*cos_dih_ideal_dev[itor]+sinphi*sin_dih_ideal_dev[itor];

    e_tor += alpha[itor] + beta[itor]*B
      + gam[itor]*(4.0*B*B*B - 3.0*B) 
      + delta[itor]*(sinphi*cos_dih_ideal_dev[itor]
		     -cosphi*sin_dih_ideal_dev[itor]);
    
  }
  
  return;
  
}

void torsional_forces()
{

  using namespace std;

  char line[2048];

  int ibead,jbead,kbead,lbead;
  double dx1, dy1, dz1;
  double dx2, dy2, dz2;
  double dx3, dy3, dz3;
  double C_11,C_12,C_13;
  double C_22,C_23;
  double C_33;
  double D_12,D_23;
  double coeff1, coeff2, coeff3;
  const static double tol = 1.0e-7;
  double cosphi,sinphi,sign,B;
  double dVdcosphi, dBdcosphi, dCdcosphi, dDdcosphi;
  coord grad[4+1]; // gradients
  double co1, co2, co3, co4, co5, co6;

  for( int itor=1; itor<=ntor; itor++ ) {
    
    ibead = ibead_tor[itor];
    jbead = jbead_tor[itor];
    kbead = kbead_tor[itor];
    lbead = lbead_tor[itor];

    dx1 = unc_pos[jbead].x - unc_pos[ibead].x;
    dy1 = unc_pos[jbead].y - unc_pos[ibead].y;
    dz1 = unc_pos[jbead].z - unc_pos[ibead].z;

    dx2 = unc_pos[kbead].x - unc_pos[jbead].x;
    dy2 = unc_pos[kbead].y - unc_pos[jbead].y;
    dz2 = unc_pos[kbead].z - unc_pos[jbead].z;

    dx3 = unc_pos[lbead].x - unc_pos[kbead].x;
    dy3 = unc_pos[lbead].y - unc_pos[kbead].y;
    dz3 = unc_pos[lbead].z - unc_pos[kbead].z;

    // min images

    dx1 -= boxl*rnd(dx1/boxl);
    dy1 -= boxl*rnd(dy1/boxl);
    dz1 -= boxl*rnd(dz1/boxl);
       
    dx2 -= boxl*rnd(dx2/boxl);
    dy2 -= boxl*rnd(dy2/boxl);
    dz2 -= boxl*rnd(dz2/boxl);
       
    dx3 -= boxl*rnd(dx3/boxl);
    dy3 -= boxl*rnd(dy3/boxl);
    dz3 -= boxl*rnd(dz3/boxl);

    C_11 = dx1*dx1 + dy1*dy1 + dz1*dz1;
    C_12 = dx1*dx2 + dy1*dy2 + dz1*dz2;
    C_13 = dx1*dx3 + dy1*dy3 + dz1*dz3;
    C_22 = dx2*dx2 + dy2*dy2 + dz2*dz2;
    C_23 = dx2*dx3 + dy2*dy3 + dz2*dz3;
    C_33 = dx3*dx3 + dy3*dy3 + dz3*dz3;

    D_12 = C_11*C_22 - C_12*C_12;
    D_23 = C_22*C_33 - C_23*C_23;    

    coeff3 = D_12*D_23;
    coeff1 = sqrt(coeff3);
    coeff2 = (C_12*C_23-C_13*C_22)/coeff1;

    cosphi = coeff2;
    sinphi = 1.0 - cosphi*cosphi;
    if( sinphi<tol ) continue; // sinphi is zero; force will blow up!!!
    sinphi = sqrt(sinphi);
    sign = (dz1*dy3-dy1*dz3)*dx2 + (dx1*dz3-dz1*dx3)*dy2
      + (dy1*dx3-dx1*dy3)*dz2;
    if( sign<0.0 ) sinphi = -sinphi;    
    
    B = cosphi*cos_dih_ideal_dev[itor] + sinphi*sin_dih_ideal_dev[itor];
    dBdcosphi = cos_dih_ideal_dev[itor] - cosphi/sinphi * sin_dih_ideal_dev[itor];
    dCdcosphi = (12*B*B - 3)*dBdcosphi;
    dDdcosphi = -cosphi/sinphi*cos_dih_ideal_dev[itor] - sin_dih_ideal_dev[itor];
    dVdcosphi = beta[itor]*dBdcosphi + gam[itor]*dCdcosphi + delta[itor]*dDdcosphi;    

    co1 = 1.0/coeff1;
    co2 = co1*co1*co1;
    co3 = co2*D_23*(C_12*C_23-C_13*C_22);
    co4 = co3*C_22;
    co5 = co1*C_23+co3*C_12;
    co6 = co1*C_22;
    
    grad[1].x = co4*dx1-co5*dx2+co6*dx3;
    grad[1].y = co4*dy1-co5*dy2+co6*dy3;
    grad[1].z = co4*dz1-co5*dz2+co6*dz3;

    co3 = co2*(C_12*C_23-C_13*C_22);
    co4 = co3*D_23*(C_12+C_22)+co1*C_23;
    co5 = co1*(C_23+2.0*C_13)+
      co3*(D_23*(C_11+C_12)+D_12*C_33);
    co6 = co1*(C_12+C_22) + co3*D_12*C_23;

    grad[2].x = -co4*dx1+co5*dx2-co6*dx3;
    grad[2].y = -co4*dy1+co5*dy2-co6*dy3;
    grad[2].z = -co4*dz1+co5*dz2-co6*dz3;

    co4 = co1*(C_22+C_23)+co3*(D_23*C_12);
    co5 = co1*(C_12+2.0*C_13)+
      co3*(D_23*C_11 + D_12*(C_33+C_23));
    co6 = co1*C_12 + co3*D_12*(C_22+C_23);

    grad[3].x = co4*dx1-co5*dx2+co6*dx3;
    grad[3].y = co4*dy1-co5*dy2+co6*dy3;
    grad[3].z = co4*dz1-co5*dz2+co6*dz3;

    co4 = co1*C_22;
    co5 = co1*C_12 + co3*D_12*C_23;
    co6 = co3*D_12*C_22;

    grad[4].x = -co4*dx1+co5*dx2-co6*dx3;
    grad[4].y = -co4*dy1+co5*dy2-co6*dy3;
    grad[4].z = -co4*dz1+co5*dz2-co6*dz3;    

    force[ibead].x -= dVdcosphi*grad[1].x;
    force[ibead].y -= dVdcosphi*grad[1].y;
    force[ibead].z -= dVdcosphi*grad[1].z;

    force[jbead].x -= dVdcosphi*grad[2].x;
    force[jbead].y -= dVdcosphi*grad[2].y;
    force[jbead].z -= dVdcosphi*grad[2].z;

    force[kbead].x -= dVdcosphi*grad[3].x;
    force[kbead].y -= dVdcosphi*grad[3].y;
    force[kbead].z -= dVdcosphi*grad[3].z;

    force[lbead].x -= dVdcosphi*grad[4].x;
    force[lbead].y -= dVdcosphi*grad[4].y;
    force[lbead].z -= dVdcosphi*grad[4].z;
    


  }

}

void vdw_energy()
{

  using namespace std;

  int ibead,jbead;
  double dx,dy,dz,d,d2,d6,d12;
  char line[2048];

  e_vdw_rr = 0.0;
  e_vdw_rr_att = 0.0;
  e_vdw_rr_rep = 0.0;

  for( int i=1; i<=nil_att; i++ ) {

    ibead = ibead_pair_list_att[i];
    jbead = jbead_pair_list_att[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);
    
    d2 = dx*dx+dy*dy+dz*dz;
    d6 = d2*d2*d2;
    d12 = d6*d6;

    e_vdw_rr_att += (pl_lj_nat_pdb_dist12[i]/d12)-2.0*(pl_lj_nat_pdb_dist6[i]/d6);

  }

  e_vdw_rr_att *= coeff_att;

  for( int i=1; i<=nil_rep; i++ ) {
    
    ibead = ibead_pair_list_rep[i];
    jbead = jbead_pair_list_rep[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;
    d6 = d2*d2*d2;
    d12 = d6*d6;    

    e_vdw_rr_rep += (sigma_rep12/d12+switch_fnb[i]*sigma_rep6/d6);

  }

  e_vdw_rr_rep *= coeff_rep;

  e_vdw_rr = e_vdw_rr_att + e_vdw_rr_rep;

  return;

}

void vdw_forces()
{

  using namespace std;

  char line[2048];

  int ibead,jbead;
  double dx,dy,dz,d,d2,d6,d12;
  double fx,fy,fz;
  double co1;
  const static double tol = 1.0e-7;
  const static double rep_tol = sigma_rep2*tol;

  for( int i=1; i<=nil_att; i++ ) {

    ibead = ibead_pair_list_att[i];
    jbead = jbead_pair_list_att[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);
    
    d2 = dx*dx+dy*dy+dz*dz;
    if( d2 < tol*lj_nat_pdb_dist2[i] ) continue;
    d6 = d2*d2*d2;
    d12 = d6*d6;
    
    co1 = force_coeff_att/d2*((pl_lj_nat_pdb_dist12[i]/d12)-(pl_lj_nat_pdb_dist6[i]/d6));

    fx = co1*dx;
    fy = co1*dy;
    fz = co1*dz;

    force[ibead].x += fx;
    force[ibead].y += fy;
    force[ibead].z += fz;

    force[jbead].x -= fx;
    force[jbead].y -= fy;
    force[jbead].z -= fz;

  }

  for( int i=1; i<=nil_rep; i++ ) {

    ibead = ibead_pair_list_rep[i];
    jbead = jbead_pair_list_rep[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;
    if( d2 <  rep_tol ) continue;    
    d6 = d2*d2*d2;
    d12 = d6*d6;    

    co1 = force_coeff_rep/d2*
      (2.0*sigma_rep12/d12+switch_fnb[i]*sigma_rep6/d6);
    
    fx = co1*dx;
    fy = co1*dy;
    fz = co1*dz;

    force[ibead].x += fx;
    force[ibead].y += fy;
    force[ibead].z += fz;

    force[jbead].x -= fx;
    force[jbead].y -= fy;
    force[jbead].z -= fz;

  }
}

void stacking_energy()
{

  using namespace std;

  // beads involded in a particular stacking interaction

  int bead1,bead2,bead3,bead4;
  int bead5,bead6,bead7,bead8;
  coord r_ang[nstck_ang+1][2+1];
  coord r_dist[nstck_dist+1];  
  coord r_tor[nstck_tor+1][3+1];
  double C_ang[nstck_ang+1][2+1][2+1]; // dot products
  double C_dist[nstck_dist+1];
  double C_tor[nstck_tor+1][3+1][3+1];
  double D_tor[nstck_tor+1][3+1][3+1]; // sums of dots
  double d_ang[nstck_ang+1][2+1]; // distances
  double d_dist[nstck_dist+1];
  double d_dist0[nstck_dist+1];
  double distM0[nstck_dist+1];
  double d_tor[nstck_tor+1][3+1];
  const int nstck_ang_plus_tor = nstck_ang + nstck_tor;
  double COS[nstck_ang_plus_tor+1];
  double SIN[nstck_ang_plus_tor+1];
  double COS0[nstck_ang_plus_tor+1];
  double SIN0[nstck_ang_plus_tor+1];
  double COSM0[nstck_ang_plus_tor+1];
  double SINM0[nstck_ang_plus_tor+1];
  double sign[nstck_tor+1]; // just used for determining the actual value of the sin

  e_stack = 0.0;
  for( int i=1; i<=nstck; i++ ) {

    bead1 = rna_stcks[i].ibead_ang[1];
    bead2 = rna_stcks[i].jbead_ang[1];
    bead3 = rna_stcks[i].kbead_ang[1];
    bead4 = rna_stcks[i].kbead_ang[2];
    bead5 = rna_stcks[i].ibead_ang[3];
    bead6 = rna_stcks[i].jbead_ang[3];
    bead7 = rna_stcks[i].kbead_ang[3];
    bead8 = rna_stcks[i].kbead_ang[4];

    // bond vectors for first angle in stack

    r_ang[1][1].x = unc_pos[bead1].x - unc_pos[bead2].x;
    r_ang[1][1].y = unc_pos[bead1].y - unc_pos[bead2].y;
    r_ang[1][1].z = unc_pos[bead1].z - unc_pos[bead2].z;
    r_ang[1][2].x = unc_pos[bead3].x - unc_pos[bead2].x;
    r_ang[1][2].y = unc_pos[bead3].y - unc_pos[bead2].y;
    r_ang[1][2].z = unc_pos[bead3].z - unc_pos[bead2].z;

    // min images

    r_ang[1][1].x -= boxl*rnd(r_ang[1][1].x/boxl);
    r_ang[1][1].y -= boxl*rnd(r_ang[1][1].y/boxl);
    r_ang[1][1].z -= boxl*rnd(r_ang[1][1].z/boxl);
    r_ang[1][2].x -= boxl*rnd(r_ang[1][2].x/boxl);
    r_ang[1][2].y -= boxl*rnd(r_ang[1][2].y/boxl);
    r_ang[1][2].z -= boxl*rnd(r_ang[1][2].z/boxl);


    // bond vectors for second angle in stack

    r_ang[2][1].x = unc_pos[bead2].x - unc_pos[bead3].x;
    r_ang[2][1].y = unc_pos[bead2].y - unc_pos[bead3].y;
    r_ang[2][1].z = unc_pos[bead2].z - unc_pos[bead3].z;
    r_ang[2][2].x = unc_pos[bead4].x - unc_pos[bead3].x;
    r_ang[2][2].y = unc_pos[bead4].y - unc_pos[bead3].y;
    r_ang[2][2].z = unc_pos[bead4].z - unc_pos[bead3].z;

    // min images

    r_ang[2][1].x -= boxl*rnd(r_ang[2][1].x/boxl);
    r_ang[2][1].y -= boxl*rnd(r_ang[2][1].y/boxl);
    r_ang[2][1].z -= boxl*rnd(r_ang[2][1].z/boxl);
    r_ang[2][2].x -= boxl*rnd(r_ang[2][2].x/boxl);
    r_ang[2][2].y -= boxl*rnd(r_ang[2][2].y/boxl);
    r_ang[2][2].z -= boxl*rnd(r_ang[2][2].z/boxl);

    // bond vectors for third angle in stack

    r_ang[3][1].x = unc_pos[bead5].x - unc_pos[bead6].x;
    r_ang[3][1].y = unc_pos[bead5].y - unc_pos[bead6].y;
    r_ang[3][1].z = unc_pos[bead5].z - unc_pos[bead6].z;
    r_ang[3][2].x = unc_pos[bead7].x - unc_pos[bead6].x;
    r_ang[3][2].y = unc_pos[bead7].y - unc_pos[bead6].y;
    r_ang[3][2].z = unc_pos[bead7].z - unc_pos[bead6].z;

    // min images

    r_ang[3][1].x -= boxl*rnd(r_ang[3][1].x/boxl);
    r_ang[3][1].y -= boxl*rnd(r_ang[3][1].y/boxl);
    r_ang[3][1].z -= boxl*rnd(r_ang[3][1].z/boxl);
    r_ang[3][2].x -= boxl*rnd(r_ang[3][2].x/boxl);
    r_ang[3][2].y -= boxl*rnd(r_ang[3][2].y/boxl);
    r_ang[3][2].z -= boxl*rnd(r_ang[3][2].z/boxl);

    // bond vectors for fourth angle in stack

    r_ang[4][1].x = unc_pos[bead6].x - unc_pos[bead7].x;
    r_ang[4][1].y = unc_pos[bead6].y - unc_pos[bead7].y;
    r_ang[4][1].z = unc_pos[bead6].z - unc_pos[bead7].z;
    r_ang[4][2].x = unc_pos[bead8].x - unc_pos[bead7].x;
    r_ang[4][2].y = unc_pos[bead8].y - unc_pos[bead7].y;
    r_ang[4][2].z = unc_pos[bead8].z - unc_pos[bead7].z;

    // min images

    r_ang[4][1].x -= boxl*rnd(r_ang[4][1].x/boxl);
    r_ang[4][1].y -= boxl*rnd(r_ang[4][1].y/boxl);
    r_ang[4][1].z -= boxl*rnd(r_ang[4][1].z/boxl);
    r_ang[4][2].x -= boxl*rnd(r_ang[4][2].x/boxl);
    r_ang[4][2].y -= boxl*rnd(r_ang[4][2].y/boxl);
    r_ang[4][2].z -= boxl*rnd(r_ang[4][2].z/boxl);

    // bond vectors for first distance in stack

    r_dist[1].x = unc_pos[bead3].x - unc_pos[bead2].x;
    r_dist[1].y = unc_pos[bead3].y - unc_pos[bead2].y;
    r_dist[1].z = unc_pos[bead3].z - unc_pos[bead2].z;

    // min images

    r_dist[1].x -= boxl*rnd(r_dist[1].x/boxl);
    r_dist[1].y -= boxl*rnd(r_dist[1].y/boxl);
    r_dist[1].z -= boxl*rnd(r_dist[1].z/boxl);

    // bond vectors for second distance in stack

    r_dist[2].x = unc_pos[bead7].x - unc_pos[bead6].x;
    r_dist[2].y = unc_pos[bead7].y - unc_pos[bead6].y;
    r_dist[2].z = unc_pos[bead7].z - unc_pos[bead6].z;

    // min images

    r_dist[2].x -= boxl*rnd(r_dist[2].x/boxl);
    r_dist[2].y -= boxl*rnd(r_dist[2].y/boxl);
    r_dist[2].z -= boxl*rnd(r_dist[2].z/boxl);

    // bond vectors for first torsion in stack

    r_tor[1][1].x = unc_pos[bead1].x - unc_pos[bead2].x;
    r_tor[1][1].y = unc_pos[bead1].y - unc_pos[bead2].y;    
    r_tor[1][1].z = unc_pos[bead1].z - unc_pos[bead2].z;    
    r_tor[1][2].x = unc_pos[bead5].x - unc_pos[bead1].x;
    r_tor[1][2].y = unc_pos[bead5].y - unc_pos[bead1].y;    
    r_tor[1][2].z = unc_pos[bead5].z - unc_pos[bead1].z;    
    r_tor[1][3].x = unc_pos[bead6].x - unc_pos[bead5].x;
    r_tor[1][3].y = unc_pos[bead6].y - unc_pos[bead5].y;    
    r_tor[1][3].z = unc_pos[bead6].z - unc_pos[bead5].z;    

    // min images

    r_tor[1][1].x -= boxl*rnd(r_tor[1][1].x/boxl);
    r_tor[1][1].y -= boxl*rnd(r_tor[1][1].y/boxl);
    r_tor[1][1].z -= boxl*rnd(r_tor[1][1].z/boxl);
    r_tor[1][2].x -= boxl*rnd(r_tor[1][2].x/boxl);
    r_tor[1][2].y -= boxl*rnd(r_tor[1][2].y/boxl);
    r_tor[1][2].z -= boxl*rnd(r_tor[1][2].z/boxl);
    r_tor[1][3].x -= boxl*rnd(r_tor[1][3].x/boxl);
    r_tor[1][3].y -= boxl*rnd(r_tor[1][3].y/boxl);
    r_tor[1][3].z -= boxl*rnd(r_tor[1][3].z/boxl);

    // bond vectors for second torsion in stack

    r_tor[2][1].x = unc_pos[bead8].x - unc_pos[bead7].x;
    r_tor[2][1].y = unc_pos[bead8].y - unc_pos[bead7].y;    
    r_tor[2][1].z = unc_pos[bead8].z - unc_pos[bead7].z;    
    r_tor[2][2].x = unc_pos[bead4].x - unc_pos[bead8].x;
    r_tor[2][2].y = unc_pos[bead4].y - unc_pos[bead8].y;    
    r_tor[2][2].z = unc_pos[bead4].z - unc_pos[bead8].z;    
    r_tor[2][3].x = unc_pos[bead3].x - unc_pos[bead4].x;
    r_tor[2][3].y = unc_pos[bead3].y - unc_pos[bead4].y;    
    r_tor[2][3].z = unc_pos[bead3].z - unc_pos[bead4].z;        

    // min images

    r_tor[2][1].x -= boxl*rnd(r_tor[2][1].x/boxl);
    r_tor[2][1].y -= boxl*rnd(r_tor[2][1].y/boxl);
    r_tor[2][1].z -= boxl*rnd(r_tor[2][1].z/boxl);
    r_tor[2][2].x -= boxl*rnd(r_tor[2][2].x/boxl);
    r_tor[2][2].y -= boxl*rnd(r_tor[2][2].y/boxl);
    r_tor[2][2].z -= boxl*rnd(r_tor[2][2].z/boxl);
    r_tor[2][3].x -= boxl*rnd(r_tor[2][3].x/boxl);
    r_tor[2][3].y -= boxl*rnd(r_tor[2][3].y/boxl);
    r_tor[2][3].z -= boxl*rnd(r_tor[2][3].z/boxl);

    // NOW STORE USEFUL DOT PRODUCTS AND SUMS OF DOT PRODUCTS

    for( int j=1; j<=nstck_ang; j++ ) {
      
      d_ang[j][1] = sqrt(r_ang[j][1].x*r_ang[j][1].x + 
	r_ang[j][1].y*r_ang[j][1].y + r_ang[j][1].z*r_ang[j][1].z);
      C_ang[j][1][2] = r_ang[j][1].x*r_ang[j][2].x + 
	r_ang[j][1].y*r_ang[j][2].y + r_ang[j][1].z*r_ang[j][2].z;
      d_ang[j][2] = sqrt(r_ang[j][2].x*r_ang[j][2].x + 
	r_ang[j][2].y*r_ang[j][2].y + r_ang[j][2].z*r_ang[j][2].z);

    }

    for( int j=1; j<=nstck_dist; j++ ) {
      
      d_dist[j] = sqrt(r_dist[j].x*r_dist[j].x +
	r_dist[j].y*r_dist[j].y + r_dist[j].z*r_dist[j].z);
      
    }

    for( int j=1; j<=nstck_tor; j++ ) {
      
      C_tor[j][1][1] = r_tor[j][1].x*r_tor[j][1].x +
	r_tor[j][1].y*r_tor[j][1].y + r_tor[j][1].z*r_tor[j][1].z;
      C_tor[j][1][2] = r_tor[j][1].x*r_tor[j][2].x +
	r_tor[j][1].y*r_tor[j][2].y + r_tor[j][1].z*r_tor[j][2].z;
      C_tor[j][1][3] = r_tor[j][1].x*r_tor[j][3].x +
	r_tor[j][1].y*r_tor[j][3].y + r_tor[j][1].z*r_tor[j][3].z;
      C_tor[j][2][2] = r_tor[j][2].x*r_tor[j][2].x +
	r_tor[j][2].y*r_tor[j][2].y + r_tor[j][2].z*r_tor[j][2].z;      
      C_tor[j][2][3] = r_tor[j][2].x*r_tor[j][3].x +
	r_tor[j][2].y*r_tor[j][3].y + r_tor[j][2].z*r_tor[j][3].z;      
      C_tor[j][3][3] = r_tor[j][3].x*r_tor[j][3].x +
	r_tor[j][3].y*r_tor[j][3].y + r_tor[j][3].z*r_tor[j][3].z;      
      d_tor[j][1] = sqrt(C_tor[j][1][1]);
      d_tor[j][2] = sqrt(C_tor[j][2][2]);
      d_tor[j][3] = sqrt(C_tor[j][3][3]);

      D_tor[j][1][2] = C_tor[j][1][1]*C_tor[j][2][2] - 
	C_tor[j][1][2]*C_tor[j][1][2];
      D_tor[j][2][3] = C_tor[j][2][2]*C_tor[j][3][3] - 
	C_tor[j][2][3]*C_tor[j][2][3];

    }

    // STORE IDEAL DIHEDRAL VALUES

    for( int j=1; j<=nstck_ang; j++ ) {
      
      COS0[j] = rna_stcks[i].cos_pdb_ang[j];
      SIN0[j] = rna_stcks[i].sin_pdb_ang[j];

    }

    for( int j=1; j<=nstck_tor; j++ ) {
    
      COS0[nstck_ang+j] = rna_stcks[i].cos_pdb_tor[j];
      SIN0[nstck_ang+j] = rna_stcks[i].sin_pdb_tor[j];
      
    }

    // STORE IDEAL DISTANCES

    for( int j=1; j<=nstck_dist; j++ ) {
      d_dist0[j] = rna_stcks[i].pdb_dist[j];
    }

    // CALCULATE ELMTS OF sign ARRAY

    for( int j=1; j<=nstck_tor; j++ ) {
    
      sign[j] = (r_tor[j][1].z*r_tor[j][3].y-r_tor[j][1].y*r_tor[j][3].z)*r_tor[j][2].x +
	(r_tor[j][1].x*r_tor[j][3].z-r_tor[j][1].z*r_tor[j][3].x)*r_tor[j][2].y +
	(r_tor[j][1].y*r_tor[j][3].x-r_tor[j][1].x*r_tor[j][3].y)*r_tor[j][2].z;

    }
    
    // CALCULATE ANGLES

    for( int j=1; j<=nstck_ang; j++ ) {
      COS[j] = C_ang[j][1][2]/d_ang[j][1]/d_ang[j][2];
      SIN[j] = sqrt(1.0 - COS[j]*COS[j]);
    }

    for( int j=1; j<=nstck_tor; j++ ) {
      COS[nstck_ang+j] = C_tor[j][1][2]*C_tor[j][2][3] - C_tor[j][1][3]*C_tor[j][2][2];
      COS[nstck_ang+j] /= sqrt(D_tor[j][1][2]*D_tor[j][2][3]);
      SIN[nstck_ang+j] = sqrt(1.0 - COS[nstck_ang+j]*COS[nstck_ang+j]);
      if( sign[j] < 0.0 ) SIN[nstck_ang+j] = -SIN[nstck_ang+j];
    }

    // CALC COSM0 AND SINM0

    for( int j=1; j<=nstck_ang_plus_tor; j++ ) {
      COSM0[j] = COS[j]*COS0[j] + SIN[j]*SIN0[j];
      SINM0[j] = SIN[j]*COS0[j] - COS[j]*SIN0[j];
    }

    // CALC distM0

    for( int j=1; j<=nstck_dist; j++ ) {
      distM0[j] = d_dist[j] - d_dist0[j];
    }

    // CALCULATE POTENTIAL
    
    e_stack += rna_stcks[i].delta_G *
      exp( -alpha_st*(SINM0[1]*SINM0[1] + SINM0[2]*SINM0[2] +
		      SINM0[3]*SINM0[3] + SINM0[4]*SINM0[4])
	   -beta_st*(distM0[1]*distM0[1] + distM0[2]*distM0[2]) 
	   -gamma_st*(SINM0[5]*SINM0[5] + SINM0[6]*SINM0[6]));

  }
  
  return;
  
}

void stacking_forces()
{

  using namespace std;

  // beads involded in a particular stacking interaction

  int bead1,bead2,bead3,bead4;
  int bead5,bead6,bead7,bead8;
  coord r_ang[nstck_ang+1][2+1];
  coord r_dist[nstck_dist+1];  
  coord r_tor[nstck_tor+1][3+1];
  double C_ang[nstck_ang+1][2+1][2+1]; // dot products
  double C_dist[nstck_dist+1];
  double C_tor[nstck_tor+1][3+1][3+1];
  double D_tor[nstck_tor+1][3+1][3+1]; // sums of dots
  double d_ang[nstck_ang+1][2+1]; // distances
  double d_dist[nstck_dist+1];
  double d_dist0[nstck_dist+1];
  double distM0[nstck_dist+1];
  double d_tor[nstck_tor+1][3+1];
  const int nstck_ang_plus_tor = nstck_ang + nstck_tor;
  double COS[nstck_ang_plus_tor+1];
  double SIN[nstck_ang_plus_tor+1];
  double SIN2[nstck_ang_plus_tor+1];
  double COS0[nstck_ang_plus_tor+1];
  double SIN0[nstck_ang_plus_tor+1];
  double COSM0[nstck_ang_plus_tor+1];
  double SINM0[nstck_ang_plus_tor+1];
  double sign[nstck_tor+1]; // just used for determining the actual value of the sin
  double pot;
  
  //FORCE DECLARATIONS
  const int nbead_stck = 8;
  const int nfunc_stck = 8;
  coord grad[nfunc_stck+1][nbead_stck+1]; // gradients
  double co1,co2,co3,co4,co5,co6;
  char line[2048];
  const static double tol = 1.0e-7;

  for( int i=1; i<=nstck; i++ ) {

    pot = 0.0;
    
    bead1 = rna_stcks[i].ibead_ang[1];
    bead2 = rna_stcks[i].jbead_ang[1];
    bead3 = rna_stcks[i].kbead_ang[1];
    bead4 = rna_stcks[i].kbead_ang[2];
    bead5 = rna_stcks[i].ibead_ang[3];
    bead6 = rna_stcks[i].jbead_ang[3];
    bead7 = rna_stcks[i].kbead_ang[3];
    bead8 = rna_stcks[i].kbead_ang[4];

    // bond vectors for first angle in stack

    r_ang[1][1].x = unc_pos[bead1].x - unc_pos[bead2].x;
    r_ang[1][1].y = unc_pos[bead1].y - unc_pos[bead2].y;
    r_ang[1][1].z = unc_pos[bead1].z - unc_pos[bead2].z;
    r_ang[1][2].x = unc_pos[bead3].x - unc_pos[bead2].x;
    r_ang[1][2].y = unc_pos[bead3].y - unc_pos[bead2].y;
    r_ang[1][2].z = unc_pos[bead3].z - unc_pos[bead2].z;

    // min images

    r_ang[1][1].x -= boxl*rnd(r_ang[1][1].x/boxl);
    r_ang[1][1].y -= boxl*rnd(r_ang[1][1].y/boxl);
    r_ang[1][1].z -= boxl*rnd(r_ang[1][1].z/boxl);
    r_ang[1][2].x -= boxl*rnd(r_ang[1][2].x/boxl);
    r_ang[1][2].y -= boxl*rnd(r_ang[1][2].y/boxl);
    r_ang[1][2].z -= boxl*rnd(r_ang[1][2].z/boxl);


    // bond vectors for second angle in stack

    r_ang[2][1].x = unc_pos[bead2].x - unc_pos[bead3].x;
    r_ang[2][1].y = unc_pos[bead2].y - unc_pos[bead3].y;
    r_ang[2][1].z = unc_pos[bead2].z - unc_pos[bead3].z;
    r_ang[2][2].x = unc_pos[bead4].x - unc_pos[bead3].x;
    r_ang[2][2].y = unc_pos[bead4].y - unc_pos[bead3].y;
    r_ang[2][2].z = unc_pos[bead4].z - unc_pos[bead3].z;

    // min images

    r_ang[2][1].x -= boxl*rnd(r_ang[2][1].x/boxl);
    r_ang[2][1].y -= boxl*rnd(r_ang[2][1].y/boxl);
    r_ang[2][1].z -= boxl*rnd(r_ang[2][1].z/boxl);
    r_ang[2][2].x -= boxl*rnd(r_ang[2][2].x/boxl);
    r_ang[2][2].y -= boxl*rnd(r_ang[2][2].y/boxl);
    r_ang[2][2].z -= boxl*rnd(r_ang[2][2].z/boxl);

    // bond vectors for third angle in stack

    r_ang[3][1].x = unc_pos[bead5].x - unc_pos[bead6].x;
    r_ang[3][1].y = unc_pos[bead5].y - unc_pos[bead6].y;
    r_ang[3][1].z = unc_pos[bead5].z - unc_pos[bead6].z;
    r_ang[3][2].x = unc_pos[bead7].x - unc_pos[bead6].x;
    r_ang[3][2].y = unc_pos[bead7].y - unc_pos[bead6].y;
    r_ang[3][2].z = unc_pos[bead7].z - unc_pos[bead6].z;

    // min images

    r_ang[3][1].x -= boxl*rnd(r_ang[3][1].x/boxl);
    r_ang[3][1].y -= boxl*rnd(r_ang[3][1].y/boxl);
    r_ang[3][1].z -= boxl*rnd(r_ang[3][1].z/boxl);
    r_ang[3][2].x -= boxl*rnd(r_ang[3][2].x/boxl);
    r_ang[3][2].y -= boxl*rnd(r_ang[3][2].y/boxl);
    r_ang[3][2].z -= boxl*rnd(r_ang[3][2].z/boxl);

    // bond vectors for fourth angle in stack

    r_ang[4][1].x = unc_pos[bead6].x - unc_pos[bead7].x;
    r_ang[4][1].y = unc_pos[bead6].y - unc_pos[bead7].y;
    r_ang[4][1].z = unc_pos[bead6].z - unc_pos[bead7].z;
    r_ang[4][2].x = unc_pos[bead8].x - unc_pos[bead7].x;
    r_ang[4][2].y = unc_pos[bead8].y - unc_pos[bead7].y;
    r_ang[4][2].z = unc_pos[bead8].z - unc_pos[bead7].z;

    // min images

    r_ang[4][1].x -= boxl*rnd(r_ang[4][1].x/boxl);
    r_ang[4][1].y -= boxl*rnd(r_ang[4][1].y/boxl);
    r_ang[4][1].z -= boxl*rnd(r_ang[4][1].z/boxl);
    r_ang[4][2].x -= boxl*rnd(r_ang[4][2].x/boxl);
    r_ang[4][2].y -= boxl*rnd(r_ang[4][2].y/boxl);
    r_ang[4][2].z -= boxl*rnd(r_ang[4][2].z/boxl);

    // bond vectors for first distance in stack

    r_dist[1].x = unc_pos[bead3].x - unc_pos[bead2].x;
    r_dist[1].y = unc_pos[bead3].y - unc_pos[bead2].y;
    r_dist[1].z = unc_pos[bead3].z - unc_pos[bead2].z;

    // min images

    r_dist[1].x -= boxl*rnd(r_dist[1].x/boxl);
    r_dist[1].y -= boxl*rnd(r_dist[1].y/boxl);
    r_dist[1].z -= boxl*rnd(r_dist[1].z/boxl);

    // bond vectors for second distance in stack

    r_dist[2].x = unc_pos[bead7].x - unc_pos[bead6].x;
    r_dist[2].y = unc_pos[bead7].y - unc_pos[bead6].y;
    r_dist[2].z = unc_pos[bead7].z - unc_pos[bead6].z;

    // min images

    r_dist[2].x -= boxl*rnd(r_dist[2].x/boxl);
    r_dist[2].y -= boxl*rnd(r_dist[2].y/boxl);
    r_dist[2].z -= boxl*rnd(r_dist[2].z/boxl);

    // bond vectors for first torsion in stack

    r_tor[1][1].x = unc_pos[bead1].x - unc_pos[bead2].x;
    r_tor[1][1].y = unc_pos[bead1].y - unc_pos[bead2].y;    
    r_tor[1][1].z = unc_pos[bead1].z - unc_pos[bead2].z;    
    r_tor[1][2].x = unc_pos[bead5].x - unc_pos[bead1].x;
    r_tor[1][2].y = unc_pos[bead5].y - unc_pos[bead1].y;    
    r_tor[1][2].z = unc_pos[bead5].z - unc_pos[bead1].z;    
    r_tor[1][3].x = unc_pos[bead6].x - unc_pos[bead5].x;
    r_tor[1][3].y = unc_pos[bead6].y - unc_pos[bead5].y;    
    r_tor[1][3].z = unc_pos[bead6].z - unc_pos[bead5].z;    

    // min images

    r_tor[1][1].x -= boxl*rnd(r_tor[1][1].x/boxl);
    r_tor[1][1].y -= boxl*rnd(r_tor[1][1].y/boxl);
    r_tor[1][1].z -= boxl*rnd(r_tor[1][1].z/boxl);
    r_tor[1][2].x -= boxl*rnd(r_tor[1][2].x/boxl);
    r_tor[1][2].y -= boxl*rnd(r_tor[1][2].y/boxl);
    r_tor[1][2].z -= boxl*rnd(r_tor[1][2].z/boxl);
    r_tor[1][3].x -= boxl*rnd(r_tor[1][3].x/boxl);
    r_tor[1][3].y -= boxl*rnd(r_tor[1][3].y/boxl);
    r_tor[1][3].z -= boxl*rnd(r_tor[1][3].z/boxl);

    // bond vectors for second torsion in stack

    r_tor[2][1].x = unc_pos[bead8].x - unc_pos[bead7].x;
    r_tor[2][1].y = unc_pos[bead8].y - unc_pos[bead7].y;    
    r_tor[2][1].z = unc_pos[bead8].z - unc_pos[bead7].z;    
    r_tor[2][2].x = unc_pos[bead4].x - unc_pos[bead8].x;
    r_tor[2][2].y = unc_pos[bead4].y - unc_pos[bead8].y;    
    r_tor[2][2].z = unc_pos[bead4].z - unc_pos[bead8].z;    
    r_tor[2][3].x = unc_pos[bead3].x - unc_pos[bead4].x;
    r_tor[2][3].y = unc_pos[bead3].y - unc_pos[bead4].y;    
    r_tor[2][3].z = unc_pos[bead3].z - unc_pos[bead4].z;        

    // min images

    r_tor[2][1].x -= boxl*rnd(r_tor[2][1].x/boxl);
    r_tor[2][1].y -= boxl*rnd(r_tor[2][1].y/boxl);
    r_tor[2][1].z -= boxl*rnd(r_tor[2][1].z/boxl);
    r_tor[2][2].x -= boxl*rnd(r_tor[2][2].x/boxl);
    r_tor[2][2].y -= boxl*rnd(r_tor[2][2].y/boxl);
    r_tor[2][2].z -= boxl*rnd(r_tor[2][2].z/boxl);
    r_tor[2][3].x -= boxl*rnd(r_tor[2][3].x/boxl);
    r_tor[2][3].y -= boxl*rnd(r_tor[2][3].y/boxl);
    r_tor[2][3].z -= boxl*rnd(r_tor[2][3].z/boxl);

    // NOW STORE USEFUL DOT PRODUCTS AND SUMS OF DOT PRODUCTS

    for( int j=1; j<=nstck_ang; j++ ) {
      
      C_ang[j][1][1] = r_ang[j][1].x*r_ang[j][1].x + 
	r_ang[j][1].y*r_ang[j][1].y + r_ang[j][1].z*r_ang[j][1].z;
      C_ang[j][1][2] = r_ang[j][1].x*r_ang[j][2].x + 
	r_ang[j][1].y*r_ang[j][2].y + r_ang[j][1].z*r_ang[j][2].z;
      C_ang[j][2][2] = r_ang[j][2].x*r_ang[j][2].x + 
	r_ang[j][2].y*r_ang[j][2].y + r_ang[j][2].z*r_ang[j][2].z;
      d_ang[j][1] = sqrt(C_ang[j][1][1]);
      d_ang[j][2] = sqrt(C_ang[j][2][2]);      

    }

    for( int j=1; j<=nstck_dist; j++ ) {
      
      C_dist[j] = r_dist[j].x*r_dist[j].x +
	r_dist[j].y*r_dist[j].y + r_dist[j].z*r_dist[j].z;
      d_dist[j] = sqrt(C_dist[j]);
      
    }

    for( int j=1; j<=nstck_tor; j++ ) {
      
      C_tor[j][1][1] = r_tor[j][1].x*r_tor[j][1].x +
	r_tor[j][1].y*r_tor[j][1].y + r_tor[j][1].z*r_tor[j][1].z;
      C_tor[j][1][2] = r_tor[j][1].x*r_tor[j][2].x +
	r_tor[j][1].y*r_tor[j][2].y + r_tor[j][1].z*r_tor[j][2].z;
      C_tor[j][1][3] = r_tor[j][1].x*r_tor[j][3].x +
	r_tor[j][1].y*r_tor[j][3].y + r_tor[j][1].z*r_tor[j][3].z;
      C_tor[j][2][2] = r_tor[j][2].x*r_tor[j][2].x +
	r_tor[j][2].y*r_tor[j][2].y + r_tor[j][2].z*r_tor[j][2].z;      
      C_tor[j][2][3] = r_tor[j][2].x*r_tor[j][3].x +
	r_tor[j][2].y*r_tor[j][3].y + r_tor[j][2].z*r_tor[j][3].z;      
      C_tor[j][3][3] = r_tor[j][3].x*r_tor[j][3].x +
	r_tor[j][3].y*r_tor[j][3].y + r_tor[j][3].z*r_tor[j][3].z;      
      d_tor[j][1] = sqrt(C_tor[j][1][1]);
      d_tor[j][2] = sqrt(C_tor[j][2][2]);
      d_tor[j][3] = sqrt(C_tor[j][3][3]);

      D_tor[j][1][2] = C_tor[j][1][1]*C_tor[j][2][2] - 
	C_tor[j][1][2]*C_tor[j][1][2];
      D_tor[j][2][3] = C_tor[j][2][2]*C_tor[j][3][3] - 
	C_tor[j][2][3]*C_tor[j][2][3];

    }

    // STORE IDEAL DIHEDRAL VALUES

    for( int j=1; j<=nstck_ang; j++ ) {
      
      COS0[j] = rna_stcks[i].cos_pdb_ang[j];
      SIN0[j] = rna_stcks[i].sin_pdb_ang[j];

    }

    for( int j=1; j<=nstck_tor; j++ ) {
    
      COS0[nstck_ang+j] = rna_stcks[i].cos_pdb_tor[j];
      SIN0[nstck_ang+j] = rna_stcks[i].sin_pdb_tor[j];
      
    }

    // STORE IDEAL DISTANCES

    for( int j=1; j<=nstck_dist; j++ ) {
      d_dist0[j] = rna_stcks[i].pdb_dist[j];
    }

    // CALCULATE ELMTS OF sign ARRAY

    for( int j=1; j<=nstck_tor; j++ ) {
    
      sign[j] = (r_tor[j][1].z*r_tor[j][3].y-r_tor[j][1].y*r_tor[j][3].z)*r_tor[j][2].x +
	(r_tor[j][1].x*r_tor[j][3].z-r_tor[j][1].z*r_tor[j][3].x)*r_tor[j][2].y +
	(r_tor[j][1].y*r_tor[j][3].x-r_tor[j][1].x*r_tor[j][3].y)*r_tor[j][2].z;

    }
    
    // CALCULATE ANGLES

    for( int j=1; j<=nstck_ang; j++ ) {
      COS[j] = C_ang[j][1][2]/d_ang[j][1]/d_ang[j][2];
      SIN2[j] = 1.0 - COS[j]*COS[j];
      if( SIN2[j]<0.0 ) SIN2[j]=0.0;
      SIN[j] = sqrt(SIN2[j]);
    }

    for( int j=1; j<=nstck_tor; j++ ) {
      COS[nstck_ang+j] = C_tor[j][1][2]*C_tor[j][2][3] - C_tor[j][1][3]*C_tor[j][2][2];
      COS[nstck_ang+j] /= sqrt(D_tor[j][1][2]*D_tor[j][2][3]);
      SIN2[nstck_ang+j] = 1.0 - COS[nstck_ang+j]*COS[nstck_ang+j];
      if( SIN2[nstck_ang+j]<0.0 ) SIN2[nstck_ang+j]=0.0;
      SIN[nstck_ang+j] = sqrt(SIN2[nstck_ang+j]);
      if( sign[j] < 0.0 ) SIN[nstck_ang+j] = -SIN[nstck_ang+j];
    }

    // CALC COSM0 AND SINM0

    for( int j=1; j<=nstck_ang_plus_tor; j++ ) {
      COSM0[j] = COS[j]*COS0[j] + SIN[j]*SIN0[j];
      SINM0[j] = SIN[j]*COS0[j] - COS[j]*SIN0[j];
    }

    // CALC distM0

    for( int j=1; j<=nstck_dist; j++ ) {
      distM0[j] = d_dist[j] - d_dist0[j];
    }

    // CALCULATE POTENTIAL
    
    pot = rna_stcks[i].delta_G *
      exp( -alpha_st*(SINM0[1]*SINM0[1] + SINM0[2]*SINM0[2] +
		      SINM0[3]*SINM0[3] + SINM0[4]*SINM0[4])
	   -beta_st*(distM0[1]*distM0[1] + distM0[2]*distM0[2]) 
	   -gamma_st*(SINM0[5]*SINM0[5] + SINM0[6]*SINM0[6]));

    // BEGIN FORCE CALCULATION

    // ZERO ALL GRADIENTS

    for( int j=1; j<=nfunc_stck; j++ ) {
      for( int k=1; k<=nbead_stck; k++ ) {
	grad[j][k].x = 0.0;
	grad[j][k].y = 0.0;
	grad[j][k].z = 0.0;
      }
    }

    // CALCULATE GRADIENTS

    if( SIN2[1]>tol ) {

      // grad of cosphi1 w.r.t. x1
      
      co1 = sqrt(C_ang[1][1][1]*C_ang[1][2][2]);
      co1 = 1.0/co1;
      co2 = co1*co1*co1*C_ang[1][1][2]*C_ang[1][2][2];
      grad[1][1].x = r_ang[1][2].x*co1 - r_ang[1][1].x*co2;
      grad[1][1].y = r_ang[1][2].y*co1 - r_ang[1][1].y*co2;
      grad[1][1].z = r_ang[1][2].z*co1 - r_ang[1][1].z*co2;
      
      
      // grad of cosphi1 w.r.t. x2
      
      co2 = co1*co1*co1;
      co3 = C_ang[1][1][2]*C_ang[1][2][2]*co2 - co1;
      co4 = C_ang[1][1][2]*C_ang[1][1][1]*co2 - co1;
      grad[1][2].x = r_ang[1][1].x*co3 + r_ang[1][2].x*co4;
      grad[1][2].y = r_ang[1][1].y*co3 + r_ang[1][2].y*co4;
      grad[1][2].z = r_ang[1][1].z*co3 + r_ang[1][2].z*co4;    
      

      // grad of cosphi1 w.r.t. x3
      
      co3 = co2*C_ang[1][1][1]*C_ang[1][1][2];
      grad[1][3].x = r_ang[1][1].x*co1 - r_ang[1][2].x*co3;
      grad[1][3].y = r_ang[1][1].y*co1 - r_ang[1][2].y*co3;
      grad[1][3].z = r_ang[1][1].z*co1 - r_ang[1][2].z*co3;

    } else {

      SIN[1] = 1.0; // dummy value

    }

    if( SIN2[2]>tol ) {

      // grad of cosphi2 w.r.t. x2
      
      co1 = sqrt(C_ang[2][1][1]*C_ang[2][2][2]);
      co1 = 1.0/co1;
      co2 = co1*co1*co1*C_ang[2][1][2]*C_ang[2][2][2];
      grad[2][2].x = r_ang[2][2].x*co1 - r_ang[2][1].x*co2;
      grad[2][2].y = r_ang[2][2].y*co1 - r_ang[2][1].y*co2;
      grad[2][2].z = r_ang[2][2].z*co1 - r_ang[2][1].z*co2;
      
      // grad of cosphi2 w.r.t. x3
      
      co2 = co1*co1*co1;
      co3 = C_ang[2][1][2]*C_ang[2][2][2]*co2 - co1;
      co4 = C_ang[2][1][2]*C_ang[2][1][1]*co2 - co1;
      grad[2][3].x = r_ang[2][1].x*co3 + r_ang[2][2].x*co4;
      grad[2][3].y = r_ang[2][1].y*co3 + r_ang[2][2].y*co4;
      grad[2][3].z = r_ang[2][1].z*co3 + r_ang[2][2].z*co4;    
      
      // grad of cosphi2 w.r.t. x4
      
      co3 = co2*C_ang[2][1][1]*C_ang[2][1][2];
      grad[2][4].x = r_ang[2][1].x*co1 - r_ang[2][2].x*co3;
      grad[2][4].y = r_ang[2][1].y*co1 - r_ang[2][2].y*co3;
      grad[2][4].z = r_ang[2][1].z*co1 - r_ang[2][2].z*co3;

    } else {

      SIN[2] = 1.0; // dummy value

    }

    if( SIN2[3]>tol ) {

      // grad of cosphi3 w.r.t. x5
      
      co1 = sqrt(C_ang[3][1][1]*C_ang[3][2][2]);
      co1 = 1.0/co1;
      co2 = co1*co1*co1*C_ang[3][1][2]*C_ang[3][2][2];
      grad[3][5].x = r_ang[3][2].x*co1 - r_ang[3][1].x*co2;
      grad[3][5].y = r_ang[3][2].y*co1 - r_ang[3][1].y*co2;
      grad[3][5].z = r_ang[3][2].z*co1 - r_ang[3][1].z*co2;
      
      // grad of cosphi3 w.r.t. x6
      
      co2 = co1*co1*co1;
      co3 = C_ang[3][1][2]*C_ang[3][2][2]*co2 - co1;
      co4 = C_ang[3][1][2]*C_ang[3][1][1]*co2 - co1;
      grad[3][6].x = r_ang[3][1].x*co3 + r_ang[3][2].x*co4;
      grad[3][6].y = r_ang[3][1].y*co3 + r_ang[3][2].y*co4;
      grad[3][6].z = r_ang[3][1].z*co3 + r_ang[3][2].z*co4;    
      
      // grad of cosphi3 w.r.t. x7
      
      co3 = co2*C_ang[3][1][1]*C_ang[3][1][2];
      grad[3][7].x = r_ang[3][1].x*co1 - r_ang[3][2].x*co3;
      grad[3][7].y = r_ang[3][1].y*co1 - r_ang[3][2].y*co3;
      grad[3][7].z = r_ang[3][1].z*co1 - r_ang[3][2].z*co3;

    } else {

      SIN[3] = 1.0; // dummy value

    }

    if( SIN2[4]>tol ) {

      // grad of cosphi4 w.r.t. x6
      
      co1 = sqrt(C_ang[4][1][1]*C_ang[4][2][2]);
      co1 = 1.0/co1;
      co2 = co1*co1*co1*C_ang[4][1][2]*C_ang[4][2][2];
      grad[4][6].x = r_ang[4][2].x*co1 - r_ang[4][1].x*co2;
      grad[4][6].y = r_ang[4][2].y*co1 - r_ang[4][1].y*co2;
      grad[4][6].z = r_ang[4][2].z*co1 - r_ang[4][1].z*co2;
      
      // grad of cosphi4 w.r.t. x7
      
      co2 = co1*co1*co1;
      co3 = C_ang[4][1][2]*C_ang[4][2][2]*co2 - co1;
      co4 = C_ang[4][1][2]*C_ang[4][1][1]*co2 - co1;
      grad[4][7].x = r_ang[4][1].x*co3 + r_ang[4][2].x*co4;
      grad[4][7].y = r_ang[4][1].y*co3 + r_ang[4][2].y*co4;
      grad[4][7].z = r_ang[4][1].z*co3 + r_ang[4][2].z*co4;    
      
      // grad of cosphi4 w.r.t. x8
      
      co3 = co2*C_ang[4][1][1]*C_ang[4][1][2];
      grad[4][8].x = r_ang[4][1].x*co1 - r_ang[4][2].x*co3;
      grad[4][8].y = r_ang[4][1].y*co1 - r_ang[4][2].y*co3;
      grad[4][8].z = r_ang[4][1].z*co1 - r_ang[4][2].z*co3;

    } else {

      SIN[4] = 1.0; // dummy value

    }
    
    // grad of r1 (func[5]) w.r.t. x3

    grad[5][3].x = r_dist[1].x/d_dist[1];
    grad[5][3].y = r_dist[1].y/d_dist[1];
    grad[5][3].z = r_dist[1].z/d_dist[1];

    // grad of r1 (func[5]) w.r.t. x2;
    
    grad[5][2].x = -grad[5][3].x;
    grad[5][2].y = -grad[5][3].y;
    grad[5][2].z = -grad[5][3].z;

    // grad of r2 (func[6]) w.r.t. x7

    grad[6][7].x = r_dist[2].x/d_dist[2];
    grad[6][7].y = r_dist[2].y/d_dist[2];
    grad[6][7].z = r_dist[2].z/d_dist[2];

    // grad of r2 (func[6]) w.r.t. x6;
    
    grad[6][6].x = -grad[6][7].x;
    grad[6][6].y = -grad[6][7].y;
    grad[6][6].z = -grad[6][7].z;

    if( SIN2[5]>tol ) {

      // grad of cospsi1 (func[7]) w.r.t x2;
      
      co1 = sqrt(D_tor[1][1][2]*D_tor[1][2][3]);
      co1 = 1.0/co1;
      co2 = co1*co1*co1;
      co3 =  co2*D_tor[1][2][3]*
	(C_tor[1][1][2]*C_tor[1][2][3]-C_tor[1][1][3]*C_tor[1][2][2]);
      co4 = co3*C_tor[1][2][2];
      co5 = co1*C_tor[1][2][3]+co3*C_tor[1][1][2];
      co6 = co1*C_tor[1][2][2];
      grad[7][2].x = co4*r_tor[1][1].x-co5*r_tor[1][2].x+co6*r_tor[1][3].x;
      grad[7][2].y = co4*r_tor[1][1].y-co5*r_tor[1][2].y+co6*r_tor[1][3].y;
      grad[7][2].z = co4*r_tor[1][1].z-co5*r_tor[1][2].z+co6*r_tor[1][3].z;
      
      // grad of cospsi1 (func[7]) w.r.t x1;
      
      co3 = co2*(C_tor[1][1][2]*C_tor[1][2][3]-C_tor[1][1][3]*C_tor[1][2][2]);
      co4 = co3*D_tor[1][2][3]*(C_tor[1][1][2]+C_tor[1][2][2])+co1*C_tor[1][2][3];
      co5 = co1*(C_tor[1][2][3]+2.0*C_tor[1][1][3])+
	co3*(D_tor[1][2][3]*(C_tor[1][1][1]+C_tor[1][1][2])+D_tor[1][1][2]*C_tor[1][3][3]);
      co6 = co1*(C_tor[1][1][2]+C_tor[1][2][2]) + co3*D_tor[1][1][2]*C_tor[1][2][3];
      grad[7][1].x = -co4*r_tor[1][1].x+co5*r_tor[1][2].x-co6*r_tor[1][3].x;
      grad[7][1].y = -co4*r_tor[1][1].y+co5*r_tor[1][2].y-co6*r_tor[1][3].y;
      grad[7][1].z = -co4*r_tor[1][1].z+co5*r_tor[1][2].z-co6*r_tor[1][3].z;
      
      // grad of cospsi1 (func[7]) w.r.t x5;
      
      co4 = co1*(C_tor[1][2][2]+C_tor[1][2][3])+co3*(D_tor[1][2][3]*C_tor[1][1][2]);
      co5 = co1*(C_tor[1][1][2]+2.0*C_tor[1][1][3])+
	co3*(D_tor[1][2][3]*C_tor[1][1][1] + D_tor[1][1][2]*(C_tor[1][3][3]+C_tor[1][2][3]));
      co6 = co1*C_tor[1][1][2] + co3*D_tor[1][1][2]*(C_tor[1][2][2]+C_tor[1][2][3]);
      grad[7][5].x = co4*r_tor[1][1].x-co5*r_tor[1][2].x+co6*r_tor[1][3].x;
      grad[7][5].y = co4*r_tor[1][1].y-co5*r_tor[1][2].y+co6*r_tor[1][3].y;
      grad[7][5].z = co4*r_tor[1][1].z-co5*r_tor[1][2].z+co6*r_tor[1][3].z;
      
      // grad of cospsi1 (func[7]) w.r.t. x6;
      
      co4 = co1*C_tor[1][2][2];
      co5 = co1*C_tor[1][1][2] + co3*D_tor[1][1][2]*C_tor[1][2][3];
      co6 = co3*D_tor[1][1][2]*C_tor[1][2][2];
      grad[7][6].x = -co4*r_tor[1][1].x+co5*r_tor[1][2].x-co6*r_tor[1][3].x;
      grad[7][6].y = -co4*r_tor[1][1].y+co5*r_tor[1][2].y-co6*r_tor[1][3].y;
      grad[7][6].z = -co4*r_tor[1][1].z+co5*r_tor[1][2].z-co6*r_tor[1][3].z;

    } else {

      SIN[5] = 1.0; // dummy value

    }
     
    if( SIN2[6]>tol ) { 

      // grad of cospsi2 (func[8]) w.r.t x7;
      
      co1 = sqrt(D_tor[2][1][2]*D_tor[2][2][3]);
      co1 = 1.0/co1;
      co2 = co1*co1*co1;
      co3 =  co2*D_tor[2][2][3]*
	(C_tor[2][1][2]*C_tor[2][2][3]-C_tor[2][1][3]*C_tor[2][2][2]);
      co4 = co3*C_tor[2][2][2];
      co5 = co1*C_tor[2][2][3]+co3*C_tor[2][1][2];
      co6 = co1*C_tor[2][2][2];
      grad[8][7].x = co4*r_tor[2][1].x-co5*r_tor[2][2].x+co6*r_tor[2][3].x;
      grad[8][7].y = co4*r_tor[2][1].y-co5*r_tor[2][2].y+co6*r_tor[2][3].y;
      grad[8][7].z = co4*r_tor[2][1].z-co5*r_tor[2][2].z+co6*r_tor[2][3].z;
      
      // grad of cospsi2 (func[8]) w.r.t x8;
      
      co3 = co2*(C_tor[2][1][2]*C_tor[2][2][3]-C_tor[2][1][3]*C_tor[2][2][2]);
      co4 = co3*D_tor[2][2][3]*(C_tor[2][1][2]+C_tor[2][2][2])+co1*C_tor[2][2][3];
      co5 = co1*(C_tor[2][2][3]+2.0*C_tor[2][1][3])+
	co3*(D_tor[2][2][3]*(C_tor[2][1][1]+C_tor[2][1][2])+D_tor[2][1][2]*C_tor[2][3][3]);
      co6 = co1*(C_tor[2][1][2]+C_tor[2][2][2]) + co3*D_tor[2][1][2]*C_tor[2][2][3];
      grad[8][8].x = -co4*r_tor[2][1].x+co5*r_tor[2][2].x-co6*r_tor[2][3].x;
      grad[8][8].y = -co4*r_tor[2][1].y+co5*r_tor[2][2].y-co6*r_tor[2][3].y;
      grad[8][8].z = -co4*r_tor[2][1].z+co5*r_tor[2][2].z-co6*r_tor[2][3].z;
      
      // grad of cospsi2 (func[8]) w.r.t x4;
      
      co4 = co1*(C_tor[2][2][2]+C_tor[2][2][3])+co3*(D_tor[2][2][3]*C_tor[2][1][2]);
      co5 = co1*(C_tor[2][1][2]+2.0*C_tor[2][1][3])+
	co3*(D_tor[2][2][3]*C_tor[2][1][1] + D_tor[2][1][2]*(C_tor[2][3][3]+C_tor[2][2][3]));
      co6 = co1*C_tor[2][1][2] + co3*D_tor[2][1][2]*(C_tor[2][2][2]+C_tor[2][2][3]);
      grad[8][4].x = co4*r_tor[2][1].x-co5*r_tor[2][2].x+co6*r_tor[2][3].x;
      grad[8][4].y = co4*r_tor[2][1].y-co5*r_tor[2][2].y+co6*r_tor[2][3].y;
      grad[8][4].z = co4*r_tor[2][1].z-co5*r_tor[2][2].z+co6*r_tor[2][3].z;
      
      // grad of cospsi2 (func[8]) w.r.t. x3;
      
      co4 = co1*C_tor[2][2][2];
      co5 = co1*C_tor[2][1][2] + co3*D_tor[2][1][2]*C_tor[2][2][3];
      co6 = co3*D_tor[2][1][2]*C_tor[2][2][2];
      grad[8][3].x = -co4*r_tor[2][1].x+co5*r_tor[2][2].x-co6*r_tor[2][3].x;
      grad[8][3].y = -co4*r_tor[2][1].y+co5*r_tor[2][2].y-co6*r_tor[2][3].y;
      grad[8][3].z = -co4*r_tor[2][1].z+co5*r_tor[2][2].z-co6*r_tor[2][3].z;

    } else {

      SIN[6] = 1.0; // dummy value

    }
      
    // CALCULATE THE FORCE ON EACH BEAD IN THE STACK

    // bead1

    force[bead1].x += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][1].x - SIN0[1]*COS[1]/SIN[1]*grad[1][1].x)
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][1].x - SIN0[5]*COS[5]/SIN[5]*grad[7][1].x);

    force[bead1].y += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][1].y - SIN0[1]*COS[1]/SIN[1]*grad[1][1].y)
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][1].y - SIN0[5]*COS[5]/SIN[5]*grad[7][1].y);

    force[bead1].z += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][1].z - SIN0[1]*COS[1]/SIN[1]*grad[1][1].z)
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][1].z - SIN0[5]*COS[5]/SIN[5]*grad[7][1].z);

    // bead2

    force[bead2].x += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][2].x - SIN0[1]*COS[1]/SIN[1]*grad[1][2].x)
      -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][2].x - SIN0[2]*COS[2]/SIN[2]*grad[2][2].x)
      -(-2.0)*beta_st*pot*distM0[1]*grad[5][2].x
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][2].x - SIN0[5]*COS[5]/SIN[5]*grad[7][2].x);

    force[bead2].y += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][2].y - SIN0[1]*COS[1]/SIN[1]*grad[1][2].y)
      -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][2].y - SIN0[2]*COS[2]/SIN[2]*grad[2][2].y)
      -(-2.0)*beta_st*pot*distM0[1]*grad[5][2].y
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][2].y - SIN0[5]*COS[5]/SIN[5]*grad[7][2].y);

    force[bead2].z += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][2].z - SIN0[1]*COS[1]/SIN[1]*grad[1][2].z)
      -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][2].z - SIN0[2]*COS[2]/SIN[2]*grad[2][2].z)
      -(-2.0)*beta_st*pot*distM0[1]*grad[5][2].z
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][2].z - SIN0[5]*COS[5]/SIN[5]*grad[7][2].z);

    // bead3

    force[bead3].x += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][3].x - SIN0[1]*COS[1]/SIN[1]*grad[1][3].x)
      -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][3].x - SIN0[2]*COS[2]/SIN[2]*grad[2][3].x)      
      -(-2.0)*beta_st*pot*distM0[1]*grad[5][3].x
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][3].x - SIN0[6]*COS[6]/SIN[6]*grad[8][3].x);      

    force[bead3].y += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][3].y - SIN0[1]*COS[1]/SIN[1]*grad[1][3].y)
      -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][3].y - SIN0[2]*COS[2]/SIN[2]*grad[2][3].y)      
      -(-2.0)*beta_st*pot*distM0[1]*grad[5][3].y
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][3].y - SIN0[6]*COS[6]/SIN[6]*grad[8][3].y);      

    force[bead3].z += -2.0*alpha_st*pot*COSM0[1]*
      (COS0[1]*grad[1][3].z - SIN0[1]*COS[1]/SIN[1]*grad[1][3].z)
      -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][3].z - SIN0[2]*COS[2]/SIN[2]*grad[2][3].z)      
      -(-2.0)*beta_st*pot*distM0[1]*grad[5][3].z
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][3].z - SIN0[6]*COS[6]/SIN[6]*grad[8][3].z);      

    // bead4

    force[bead4].x += -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][4].x - SIN0[2]*COS[2]/SIN[2]*grad[2][4].x)      
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][4].x - SIN0[6]*COS[6]/SIN[6]*grad[8][4].x);      


    force[bead4].y += -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][4].y - SIN0[2]*COS[2]/SIN[2]*grad[2][4].y)      
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][4].y - SIN0[6]*COS[6]/SIN[6]*grad[8][4].y);      

    force[bead4].z += -2.0*alpha_st*pot*COSM0[2]*
      (COS0[2]*grad[2][4].z - SIN0[2]*COS[2]/SIN[2]*grad[2][4].z)      
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][4].z - SIN0[6]*COS[6]/SIN[6]*grad[8][4].z);      

    // bead5

    force[bead5].x += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][5].x - SIN0[3]*COS[3]/SIN[3]*grad[3][5].x)
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][5].x - SIN0[5]*COS[5]/SIN[5]*grad[7][5].x);

    force[bead5].y += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][5].y - SIN0[3]*COS[3]/SIN[3]*grad[3][5].y)
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][5].y - SIN0[5]*COS[5]/SIN[5]*grad[7][5].y);

    force[bead5].z += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][5].z - SIN0[3]*COS[3]/SIN[3]*grad[3][5].z)
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][5].z - SIN0[5]*COS[5]/SIN[5]*grad[7][5].z);

    // bead6

    force[bead6].x += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][6].x - SIN0[3]*COS[3]/SIN[3]*grad[3][6].x)
      -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][6].x - SIN0[4]*COS[4]/SIN[4]*grad[4][6].x)
      -(-2.0)*beta_st*pot*distM0[2]*grad[6][6].x
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][6].x - SIN0[5]*COS[5]/SIN[5]*grad[7][6].x);

    force[bead6].y += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][6].y - SIN0[3]*COS[3]/SIN[3]*grad[3][6].y)
      -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][6].y - SIN0[4]*COS[4]/SIN[4]*grad[4][6].y)
      -(-2.0)*beta_st*pot*distM0[2]*grad[6][6].y
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][6].y - SIN0[5]*COS[5]/SIN[5]*grad[7][6].y);

    force[bead6].z += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][6].z - SIN0[3]*COS[3]/SIN[3]*grad[3][6].z)
      -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][6].z - SIN0[4]*COS[4]/SIN[4]*grad[4][6].z)
      -(-2.0)*beta_st*pot*distM0[2]*grad[6][6].z
      -2.0*gamma_st*pot*COSM0[5]*
      (COS0[5]*grad[7][6].z - SIN0[5]*COS[5]/SIN[5]*grad[7][6].z);

    // bead7

    force[bead7].x += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][7].x - SIN0[3]*COS[3]/SIN[3]*grad[3][7].x)
      -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][7].x - SIN0[4]*COS[4]/SIN[4]*grad[4][7].x)
      -(-2.0)*beta_st*pot*distM0[2]*grad[6][7].x
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][7].x - SIN0[6]*COS[6]/SIN[6]*grad[8][7].x);

    force[bead7].y += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][7].y - SIN0[3]*COS[3]/SIN[3]*grad[3][7].y)
      -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][7].y - SIN0[4]*COS[4]/SIN[4]*grad[4][7].y)
      -(-2.0)*beta_st*pot*distM0[2]*grad[6][7].y
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][7].y - SIN0[6]*COS[6]/SIN[6]*grad[8][7].y);

    force[bead7].z += -2.0*alpha_st*pot*COSM0[3]*
      (COS0[3]*grad[3][7].z - SIN0[3]*COS[3]/SIN[3]*grad[3][7].z)
      -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][7].z - SIN0[4]*COS[4]/SIN[4]*grad[4][7].z)
      -(-2.0)*beta_st*pot*distM0[2]*grad[6][7].z
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][7].z - SIN0[6]*COS[6]/SIN[6]*grad[8][7].z);

    // bead8
    
    force[bead8].x += -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][8].x - SIN0[4]*COS[4]/SIN[4]*grad[4][8].x)
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][8].x - SIN0[6]*COS[6]/SIN[6]*grad[8][8].x);      

    force[bead8].y += -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][8].y - SIN0[4]*COS[4]/SIN[4]*grad[4][8].y)
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][8].y - SIN0[6]*COS[6]/SIN[6]*grad[8][8].y);      

    force[bead8].z += -2.0*alpha_st*pot*COSM0[4]*
      (COS0[4]*grad[4][8].z - SIN0[4]*COS[4]/SIN[4]*grad[4][8].z)
      -2.0*gamma_st*pot*COSM0[6]*
      (COS0[6]*grad[8][8].z - SIN0[6]*COS[6]/SIN[6]*grad[8][8].z);      
    
  }
  
}

void electrostatic_energy()
{

  using namespace std;

  int ibead,jbead;
  double dx,dy,dz,d2,d;

  e_elec = 0.0;
  for( int i=1; i<=nelec; i++ ) {

    ibead = ibead_elec[i];
    jbead = jbead_elec[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;
    d = sqrt(d2);

    e_elec += 1.0/d*exp(-d/debye);
    
  }
  
  e_elec *= eelec_coeff;

  return;

}

void electrostatic_forces()
{

  using namespace std;

  char line[2048];

  int ibead, jbead;
  double dx,dy,dz,d2,d;
  double fx,fy,fz;
  double fact;

  for( int i=1; i<=nelec; i++ ) {

    ibead = ibead_elec[i];
    jbead = jbead_elec[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;
    d = sqrt(d2);

    fact = felec_coeff*(1.0/d2)*exp(-d/debye)*(1.0/debye+1.0/d);

    fx = fact*dx;
    fy = fact*dy;
    fz = fact*dz;

    force[ibead].x += fx;
    force[ibead].y += fy;
    force[ibead].z += fz;

    force[jbead].x -= fx;
    force[jbead].y -= fy;
    force[jbead].z -= fz;

    
  }
  
}


