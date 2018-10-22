#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_legendre.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include "Coefficients.h"  //WGS84 gravity model
#define _u_ 3.986012e14
#define _a_e_ 6378351.0
#define _a_ orbit6elem[0]
#define _e_ orbit6elem[1]
#define _i_ orbit6elem[2]
#define _O_ orbit6elem[3]
#define _w_ orbit6elem[4]
#define _f_ orbit6elem[5]
#define _Vx StVector[0]
#define _Vy StVector[1]
#define _Vz StVector[2]
#define _Rx StVector[3]
#define _Ry StVector[4]
#define _Rz StVector[5]
#define sidereal_time_of_greenwich_at_T0 0.0
#define Defination 3  //Defination must be lower than 71

void oblateness_perturbation(gsl_vector * obper, const gsl_vector * Position,double t){
  double r = gsl_blas_dnrm2(Position);
  double p = gsl_hypot(gsl_vector_get(Position,0),gsl_vector_get(Position,1));
  double x = gsl_vector_get(Position,2);
  gsl_vector * tmp = gsl_vector_alloc(3);gsl_vector_set_zero(tmp);
  gsl_vector * cache = gsl_vector_alloc(3);gsl_vector_set_zero(cache);
  double Lambda = atan2(gsl_vector_get(Position,1),gsl_vector_get(Position,0)) - (7.292e-5 * t + sidereal_time_of_greenwich_at_T0);
  double SinPhi = gsl_vector_get(Position,2)/r;
  double Sum = 0,Coe = 0;
  gsl_vector_set(tmp,2,1);
  gsl_vector_memcpy(obper,Position);
  for (int n=2;n<Defination;++n){
    Sum = 0;
    for (int m=0;m<n;++m){
      Sum -= ( (n+m+1)*gsl_sf_legendre_Plm(n,m,SinPhi) + x/p*gsl_sf_legendre_Plm(n,m+1,SinPhi) )
	     * (C[n][m]*cos(m*Lambda) + S[n][m]*sin(m*Lambda));
    }
    Sum -= (n+n+1)*gsl_sf_legendre_Plm(n,n,SinPhi)
           * (C[n][n]*cos(n*Lambda) + S[n][n]*sin(n*Lambda));
    Coe += Sum * gsl_pow_int(_a_e_/r,n);
  }
  Coe *= _u_/gsl_pow_3(r);
  gsl_vector_scale(obper,Coe);
  Sum = 0,Coe = 0;
  for (int n=2;n<Defination;++n){
    for (int m=0;m<n;++m){
      Sum += r/p*gsl_sf_legendre_Plm(n,m+1,SinPhi)
	* (C[n][m]*cos(m*Lambda) + S[n][m]*sin(m*Lambda));
    }
    Coe += Sum * gsl_pow_int(_a_e_/r,n);
  }
  Coe *= _u_/gsl_pow_2(r);
  gsl_vector_scale(tmp,Coe);
  gsl_vector_add(obper,tmp);
  for (int n=2;n<Defination;++n){
    for (int m=1;m<=n;++m){
      gsl_vector_set_zero(tmp);
      x = m * (C[n][m]*cos((m-1)*Lambda) + S[n][m]*sin((m-1)*Lambda));
      gsl_vector_set(tmp,0,x);
      x = -m * (C[n][m]*sin((m-1)*Lambda) - S[n][m]*cos((m-1)*Lambda));
      gsl_vector_set(tmp,1,x);
      gsl_vector_scale(tmp,r/p*gsl_sf_legendre_Plm(n,m,SinPhi)*gsl_pow_int(_a_e_/r,n));
      gsl_vector_add(cache,tmp);
    }
  }
  gsl_vector_scale(cache,_u_/gsl_pow_2(r));
  gsl_vector_add(obper,cache);
  gsl_vector_free(tmp);gsl_vector_free(cache);
}
void atmosphere_perturbation(gsl_vector * atmper, const gsl_vector * Velocity){
  double Rho = 1.534e-12;        //大气密度
  double C_D = 2.4;              //阻力系数
  double A = 9;                  //航天器迎风面积
  double m = 10000;              //航天器质量
  double v = gsl_blas_dnrm2(Velocity);
  gsl_vector_memcpy(atmper,Velocity);
  gsl_vector_scale(atmper,-0.5 * C_D * A / m * Rho * v);
}
int erti(double t, const double y[], double f[], void * params){
  (void)t;
  (void)params;
  gsl_vector_const_view sta = gsl_vector_const_view_array(y,6);
  gsl_vector_const_view Velocity = gsl_vector_const_subvector(&sta.vector,0,3);
  gsl_vector_const_view Position = gsl_vector_const_subvector(&sta.vector,3,3);
  gsl_vector * Acceleration = gsl_vector_alloc(3);
  gsl_vector * obper = gsl_vector_alloc(3);oblateness_perturbation(obper,&Position.vector,t);
  gsl_vector * atmper = gsl_vector_alloc(3);atmosphere_perturbation(atmper,&Velocity.vector);
  gsl_vector_memcpy(Acceleration,&Position.vector);
  gsl_vector_scale(Acceleration,-_u_/gsl_pow_3(gsl_blas_dnrm2(&Position.vector)));
  gsl_vector_add(Acceleration,obper);
  gsl_vector_add(Acceleration,atmper);
  for(int i = 0;i<3;++i) f[i] = gsl_vector_get(Acceleration,i);
  for(int i = 3;i<6;++i) f[i] = gsl_vector_get(&Velocity.vector,i-3);
  gsl_vector_free(Acceleration);gsl_vector_free(obper);gsl_vector_free(atmper);
  return GSL_SUCCESS;
}
void StVector2elem(double * orbit6elem, double * StVector){
  (void)orbit6elem;
  (void)StVector;
  gsl_vector_view sta = gsl_vector_view_array(StVector,6);
  gsl_vector_view Velocity = gsl_vector_subvector(&sta.vector,0,3);
  gsl_vector_view Position = gsl_vector_subvector(&sta.vector,3,3);
  gsl_vector * tmp1 = gsl_vector_alloc(3);gsl_vector_memcpy(tmp1,&Velocity.vector);
  gsl_vector * tmp2 = gsl_vector_alloc(3);gsl_vector_memcpy(tmp2,&Position.vector);
  double v = gsl_blas_dnrm2(&Velocity.vector);
  double r = gsl_blas_dnrm2(&Position.vector);
  double Hx = _Ry*_Vz-_Rz*_Vy,Hy = _Rz*_Vx-_Rx*_Vz,Hz = _Rx*_Vy-_Ry*_Vx;
  double h = gsl_hypot3(Hx,Hy,Hz);
  double Nx = -Hy, Ny = Hx, N = gsl_hypot(Nx,Ny);
  double rVr = 0, Ne = 0, er = 0;
  gsl_blas_ddot(&Velocity.vector,&Position.vector,&rVr);
  gsl_vector_scale(tmp1,rVr/_u_);gsl_vector_scale(tmp2,gsl_pow_2(v)/_u_-1/r);gsl_vector_sub(tmp2,tmp1);
  gsl_vector_set(tmp1,0,Nx);gsl_vector_set(tmp1,1,Ny);gsl_vector_set(tmp1,2,0);
  gsl_blas_ddot(tmp2,tmp1,&Ne);
  gsl_blas_ddot(tmp2,&Position.vector,&er);
  _i_ = acos(Hz/h);
  _O_ = (Ny >= 0) ? acos(Nx/N) : 2*M_PI-acos(Nx/N);
  _e_ = gsl_blas_dnrm2(tmp2);
  _w_ = (gsl_vector_get(tmp2,2)>=0) ? acos(Ne/N/_e_) : 2*M_PI-acos(Ne/N/_e_);
  _f_ = (rVr >= 0) ? acos(er/_e_/r) : 2*M_PI-acos(er/_e_/r);
  _a_ = gsl_pow_2(h)/(_u_-_u_*gsl_pow_2(_e_));
  gsl_vector_free(tmp1);gsl_vector_free(tmp2);
}
void elem2stVector(double * StVector, double * orbit6elem){
  (void)StVector;
  (void)orbit6elem;
  double rx[] = {cos(_f_),sin(_f_),0};
  double vx[] = {-sin(_f_),_e_+cos(_f_),0};
  double p = _a_*(1-gsl_pow_2(_e_));
  gsl_vector_view Rx_bar = gsl_vector_view_array(rx,3);
  gsl_vector_view Vx_bar = gsl_vector_view_array(vx,3);
  gsl_vector_scale(&Rx_bar.vector,p/(1+_e_*cos(_f_)));
  gsl_vector_scale(&Vx_bar.vector,sqrt(_u_/p));
  double Q_x_bar2x[] = {cos(_O_)*cos(_w_)-sin(_O_)*sin(_w_)*cos(_i_), -cos(_O_)*sin(_w_)-sin(_O_)*cos(_w_)*cos(_i_), sin(_O_)*sin(_i_),
                        sin(_O_)*cos(_w_)+cos(_O_)*sin(_w_)*cos(_i_), -sin(_O_)*sin(_w_)+cos(_O_)*cos(_w_)*cos(_i_), -cos(_O_)*sin(_i_),
                        sin(_i_)*sin(_w_),                            sin(_i_)*cos(_w_),                             cos(_i_)};
  gsl_matrix_view Qxx = gsl_matrix_view_array(Q_x_bar2x,3,3);
  gsl_vector * Rx = gsl_vector_alloc(3);
  gsl_vector * Vx = gsl_vector_alloc(3);
  gsl_blas_dgemv(CblasNoTrans,1.0,&Qxx.matrix,&Rx_bar.vector,0,Rx);
  gsl_blas_dgemv(CblasNoTrans,1.0,&Qxx.matrix,&Vx_bar.vector,0,Vx);
  _Vx = gsl_vector_get(Vx,0);_Vy = gsl_vector_get(Vx,1);_Vz = gsl_vector_get(Vx,2);
  _Rx = gsl_vector_get(Rx,0);_Ry = gsl_vector_get(Rx,1);_Rz = gsl_vector_get(Rx,2);
  gsl_vector_free(Rx);gsl_vector_free(Vx);
}
int main(){
  gsl_odeiv2_system sys = {erti,NULL,6,NULL};
  gsl_odeiv2_driver * dri = gsl_odeiv2_driver_alloc_y_new(&sys,
							  gsl_odeiv2_step_rk8pd,
							  1e-3,1e-12,1e-8);
  double t = 0.0;
  double orbit6elem[] = {6781000, 0.001, 51/57.3, 0, 0, 0};
  double StVector[6];
  //double T = 2 * M_PI * _a_ * sqrt(_a_/_u_);
  FILE * Orbit = fopen("./Orbit6elem.txt","a+");
  FILE * StVec = fopen("./StVector.txt","a+");
  elem2stVector(StVector,orbit6elem);
  int S;
  fprintf(StVec,"%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
	   t,_Vx,_Vy,_Vz,_Rx,_Ry,_Rz,
	   gsl_pow_2(gsl_hypot3(_Vx,_Vy,_Vz))/2-_u_/gsl_hypot3(_Rx,_Ry,_Rz),
	   _Ry*_Vz-_Rz*_Vy,_Rz*_Vx-_Rx*_Vz,_Rx*_Vy-_Ry*_Vx);
  fprintf(Orbit,"%.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
	    t,_a_,_e_,_i_,_O_,_w_,_f_);
  while (t < 3600*24 && gsl_hypot3(_Rx,_Ry,_Rz) > 6.371004e6) {
    S = gsl_odeiv2_driver_apply_fixed_step(dri, &t, 1, 10, StVector);
    fprintf(StVec,"%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
	   t,_Vx,_Vy,_Vz,_Rx,_Ry,_Rz,
	   gsl_pow_2(gsl_hypot3(_Vx,_Vy,_Vz))/2-_u_/gsl_hypot3(_Rx,_Ry,_Rz),
	   _Ry*_Vz-_Rz*_Vy,_Rz*_Vx-_Rx*_Vz,_Rx*_Vy-_Ry*_Vx);
    StVector2elem(orbit6elem,StVector);
    fprintf(Orbit,"%.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
	    t,_a_,_e_,_i_,_O_,_w_,_f_);
    if (S != GSL_SUCCESS){
      printf("err:driver returned %d\n",S);
      break;
    }
  }
  gsl_odeiv2_driver_free(dri);
  fclose(Orbit);Orbit=NULL;
  fclose(StVec);StVec=NULL;
  return S;
}
