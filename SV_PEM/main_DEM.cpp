//おまじない
//(意味: ここから下は普通のC++の構文で書いてありますよ)
#pragma unmanaged

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "PEMstruct.h"

#define DIM 2 //次元
#define PARTICLE_DISTANCE    0.04 //初期粒子間距離
#define HEIGHT 1.0
#define WIDTH 1.0

//type
#define GHOST  -1
#define FLUID   0
#define WALL    2
#define DUMMY_WALL  3
#define GHOST_OR_DUMMY  -1

#define SURFACE_PARTICLE 1    
#define INNER_PARTICLE   0   

#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0 
#define DIRICHLET_BOUNDARY_IS_CONNECTED     1 
#define DIRICHLET_BOUNDARY_IS_CHECKED       2   

#define ON              1
#define OFF             0

#define ARRAY_SIZE 5000	//最大粒子数
#define COEFFICIENT_OF_RESTITUTION 0.2 //剛体衝突の反発係数
#define THRESHOLD_RATIO_OF_NUMBER_DENSITY  0.97   
#define RELAXATION_COEFFICIENT_FOR_PRESSURE 0.2	//圧力計算の際の緩和係数
#define COMPRESSIBILITY (0.45E-9) //流体の圧縮率
#define RADIUS_FOR_NUMBER_DENSITY  (2.1*PARTICLE_DISTANCE) 
#define RADIUS_FOR_GRADIENT        (2.1*PARTICLE_DISTANCE) 
#define RADIUS_FOR_LAPLACIAN       (3.1*PARTICLE_DISTANCE) 
#define COLLISION_DISTANCE         (0.5*PARTICLE_DISTANCE)
#define FLUID_DENSITY        1000.0 
#define EPS             (0.01 * PARTICLE_DISTANCE)   


const double PI = 3.1415926535;		//円周率
//const double G = 0;				    //重力加速度
const double G = 9.80665;		    //重力加速度
const double rho = 10;		        //粒子の密度
double kinematic_viscosity = 0.000001; //動粘性係数

int    NumberOfParticles;
double Re_forNumberDensity,Re2_forNumberDensity; 
double Re_forGradient,     Re2_forGradient; 
double Re_forLaplacian,    Re2_forLaplacian; 
double N0_forNumberDensity;
double N0_forGradient;
double N0_forLaplacian;
double Lambda;
double collisionDistance,collisionDistance2;
double FluidDensity;

void initializeParticlePositionAndVelocity_for2dim(PEM *pe);
void calConstantParameter(void);
void calNZeroAndLambda(void);

double weight( double distance, double re );

void calGravity(PEM *pe);
void calViscosity(PEM *pe);
void moveParticle(PEM *pe, double dt);
void collision(PEM *pe, double dt);
void calPressure(PEM *pe, double dt);
void calPressureGradient(PEM *pe);
void moveParticleUsingPressureGradient(PEM *pe, double dt);

void calNumberDensity(PEM *pe);
void setBoundaryCondition(PEM *pe);
void setSourceTerm(PEM *pe, double dt);
void setMatrix(PEM *pe, double dt);
void solveSimultaniousEquationsByGaussEliminationMethod(PEM *pe);
void removeNegativePressure(PEM *pe);
void setMinimumPressure(PEM *pe);

void checkBoundaryCondition(PEM *pe);
void increaseDiagonalTerm(PEM *pe);
void exceptionalProcessingForBoundaryCondition(PEM *pe);





//系の中の各要素の情報を細かく設定する関数
//この関数を呼ぶと粒子とかが初期状態になる
void Config(int rows,int pe_in_row, int wa_n, PEM *pe, double kv){
	
	int pe_n = rows * pe_in_row;
	int n=pe_n+wa_n;
	double radius = (double)20/pe_in_row;

	//kinematic_viscosity = kv;

	initializeParticlePositionAndVelocity_for2dim(pe);
	calConstantParameter();
}

//1ステップの動きを計算する関数
void CalcStep(int pe_n, int wa_n, PEM *pe, double dt, double lc, double *param){

	calGravity(pe);
    calViscosity(pe);
    moveParticle(pe, dt);
    collision(pe, dt);
    calPressure(pe, dt);
    calPressureGradient(pe);
    moveParticleUsingPressureGradient(pe, dt);
}


void initializeParticlePositionAndVelocity_for2dim(PEM *pe){
  int iX, iY;
  int nX, nY;
  double x, y;
  int i = 0;
  int flagOfParticleGeneration;

  nX = (int)(WIDTH/PARTICLE_DISTANCE)+5;  
  nY = (int)(HEIGHT/PARTICLE_DISTANCE)+5;
  #pragma omp parallel for  num_threads(8)
  for(iX= -4;iX<nX;iX++){
    for(iY= -4;iY<nY;iY++){
      x = PARTICLE_DISTANCE * (double)(iX);
      y = PARTICLE_DISTANCE * (double)(iY);
      flagOfParticleGeneration = OFF;

      /* dummy wall region */
      if( ((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=WIDTH+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-4.0*PARTICLE_DISTANCE+EPS )&&(y<=HEIGHT+EPS)) ){  
	pe[i].type=DUMMY_WALL;
	flagOfParticleGeneration = ON;
      }

      /* wall region */
      if( ((x>-2.0*PARTICLE_DISTANCE+EPS)&&(x<=WIDTH+2.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-2.0*PARTICLE_DISTANCE+EPS )&&(y<=HEIGHT+EPS)) ){ 
	pe[i].type=WALL;
	flagOfParticleGeneration = ON;
      }

      /* wall region */
      if( ((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=WIDTH+4.0*PARTICLE_DISTANCE+EPS))&&( (y>HEIGHT-2.0*PARTICLE_DISTANCE+EPS )&&(y<=HEIGHT+EPS)) ){ 
	pe[i].type=WALL;
	flagOfParticleGeneration = ON;
      }

      /* empty region */
      if( ((x>0.0+EPS)&&(x<=WIDTH+EPS))&&( y>0.0+EPS )){  
	flagOfParticleGeneration = OFF;
      }

      /* fluid region */
      if( ((x>0.0+EPS)&&(x<=WIDTH*0.2+EPS)) &&((y>0.0+EPS)&&(y<=HEIGHT*0.8+EPS)) ){  
	pe[i].type=FLUID;
	flagOfParticleGeneration = ON;
      }

      if( flagOfParticleGeneration == ON){
	pe[i].x=x; pe[i].y=y;
	i++;
      }
    }
  }
  NumberOfParticles = i;
  for(i=0;i<NumberOfParticles;i++) {
	  pe[i].x = pe[i].x - WIDTH/2.0 - PARTICLE_DISTANCE/2.0;
	pe[i].y = pe[i].y - HEIGHT/2.0 - PARTICLE_DISTANCE/2.0;
	pe[i].vx = 0.0;
	pe[i].vy = 0.0;
	pe[i].r = PARTICLE_DISTANCE / 2.0;
  }
}


void calConstantParameter( void ){

  Re_forNumberDensity  = RADIUS_FOR_NUMBER_DENSITY;  
  Re_forGradient       = RADIUS_FOR_GRADIENT;  
  Re_forLaplacian      = RADIUS_FOR_LAPLACIAN;  
  Re2_forNumberDensity = Re_forNumberDensity*Re_forNumberDensity;
  Re2_forGradient      = Re_forGradient*Re_forGradient;
  Re2_forLaplacian     = Re_forLaplacian*Re_forLaplacian;
  calNZeroAndLambda();
  FluidDensity       = FLUID_DENSITY;
  collisionDistance  = COLLISION_DISTANCE; 
  collisionDistance2 = collisionDistance*collisionDistance;
}


void calNZeroAndLambda( void ){
  int iX, iY;
  double xj, yj, distance, distance2;
  double xi, yi;

  N0_forNumberDensity = 0.0;
  N0_forGradient      = 0.0;
  N0_forLaplacian     = 0.0;
  Lambda              = 0.0;
  xi = 0.0;  yi = 0.0;

  for(iX= -4;iX<5;iX++){
    for(iY= -4;iY<5;iY++){
	if((iX==0)&&(iY==0))continue;
	xj = PARTICLE_DISTANCE * (double)(iX);
	yj = PARTICLE_DISTANCE * (double)(iY);
	distance2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi);
	distance = sqrt(distance2);
	N0_forNumberDensity += weight(distance, Re_forNumberDensity);
	N0_forGradient      += weight(distance, Re_forGradient);
	N0_forLaplacian     += weight(distance, Re_forLaplacian);
	Lambda              += distance2 * weight(distance, Re_forLaplacian);
    }
  }
  Lambda = Lambda/N0_forLaplacian;
}


double weight( double distance, double re ){
  double weightIJ;

  if( distance >= re ){
    weightIJ = 0.0;
  }else{
    weightIJ = (re/distance) - 1.0;
  }
  return weightIJ;
}


void calGravity(PEM *pe){
  int i;

  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].type == FLUID){
	  pe[i].ax = 0.0;
      pe[i].ay = -G;
    }else{
	  pe[i].ax = 0.0;
	  pe[i].ax = 0.0;
    }
  }
}


void calViscosity(PEM *pe){
  int i,j;
  double viscosityTerm_x, viscosityTerm_y;
  double distance, distance2;
  double w;
  double xij, yij;
  double a;

  a = (kinematic_viscosity)*(2.0*DIM)/(N0_forLaplacian*Lambda);
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].type != FLUID) continue;
    viscosityTerm_x = 0.0;  viscosityTerm_y = 0.0;

    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (pe[j].type==GHOST) ) continue;
      xij = pe[j].x - pe[i].x;
      yij = pe[j].y - pe[i].y;
      distance2 = (xij*xij) + (yij*yij);
      distance = sqrt(distance2);
      if(distance<Re_forLaplacian){
	w =  weight(distance, Re_forLaplacian);
	viscosityTerm_x +=(pe[j].vx-pe[i].vx)*w;
	viscosityTerm_y +=(pe[j].vy-pe[i].vy)*w;
      }
    }
    viscosityTerm_x = viscosityTerm_x * a;
    viscosityTerm_y = viscosityTerm_y * a;
    pe[i].ax += viscosityTerm_x;
    pe[i].ay += viscosityTerm_y;
  }
}


void moveParticle(PEM *pe, double dt){
  int i;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].type == FLUID){
      pe[i].vx += pe[i].ax*dt; 
      pe[i].vy += pe[i].ay*dt; 

      pe[i].x += pe[i].vx*dt; 
      pe[i].y += pe[i].vy*dt; 
    }
    pe[i].ax=0.0;
    pe[i].ay=0.0;
  }
}


void collision(PEM *pe, double dt){
  int    i,j;
  double xij, yij;
  double distance,distance2;
  double forceDT; /* forceDT is the impulse of collision between particles */
  double mi, mj;
  double velocity_ix, velocity_iy;
  double e = COEFFICIENT_OF_RESTITUTION;
  static double VelocityAfterCollision[DIM*ARRAY_SIZE];
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){ 
    VelocityAfterCollision[i] = pe[i].vx;
    VelocityAfterCollision[i+1] = pe[i].vy;
  }
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].type == FLUID){
      mi = FluidDensity;
      velocity_ix = pe[i].vx;  
      velocity_iy = pe[i].vy;
      for(j=0;j<NumberOfParticles;j++){
		if( (j==i) || (pe[j].type==GHOST) ) continue;
		xij = pe[j].x - pe[i].x;
		yij = pe[j].y - pe[i].y;
		distance2 = (xij*xij) + (yij*yij);
		if(distance2<collisionDistance2){
		  distance = sqrt(distance2);
		  forceDT = (velocity_ix-pe[j].vx)*(xij/distance)
				   +(velocity_iy-pe[j].vy)*(yij/distance);
		  if(forceDT > 0.0){
			mj = FluidDensity;
			forceDT *= (1.0+e)*mi*mj/(mi+mj);
			velocity_ix -= (forceDT/mi)*(xij/distance); 
			velocity_iy -= (forceDT/mi)*(yij/distance);
			/*
			if(j>i){ fprintf(stderr,"WARNING: Collision occured between %d and %d particles.\n",i,j); }
			*/
		  }
		}
      }
      VelocityAfterCollision[i*DIM  ] = velocity_ix; 
      VelocityAfterCollision[i*DIM+1] = velocity_iy;
    }
  }
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].type == FLUID){
      pe[i].x += (VelocityAfterCollision[i*DIM  ]-pe[i].vx)*dt; 
      pe[i].y += (VelocityAfterCollision[i*DIM+1]-pe[i].vy)*dt;
      pe[i].vx = VelocityAfterCollision[i*DIM  ]; 
      pe[i].vy = VelocityAfterCollision[i*DIM+1];
    }
  }
}


void calPressure(PEM *pe, double dt){
  calNumberDensity(pe);
  setBoundaryCondition(pe);
  setSourceTerm(pe, dt);
  setMatrix(pe, dt);
  solveSimultaniousEquationsByGaussEliminationMethod(pe);
  removeNegativePressure(pe);
  setMinimumPressure(pe);
}


void calNumberDensity(PEM *pe){
  int    i,j;
  double xij, yij;
  double distance, distance2;
  double w;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    pe[i].numdensity = 0.0;
    if(pe[i].type == GHOST) continue;
    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (pe[j].type==GHOST) ) continue;
      xij = pe[j].x - pe[i].x;
      yij = pe[j].y - pe[i].y;
      distance2 = (xij*xij) + (yij*yij);
      distance = sqrt(distance2);
      w =  weight(distance, Re_forNumberDensity);
      pe[i].numdensity += w;
    }
  }
}


void setBoundaryCondition(PEM *pe){
  int i;
  double n0 = N0_forNumberDensity;
  double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].type==GHOST || pe[i].type== DUMMY_WALL ){
      pe[i].boundary=GHOST_OR_DUMMY;
    }else if( pe[i].numdensity < beta * n0 ){
      pe[i].boundary=SURFACE_PARTICLE;
    }else{
      pe[i].boundary=INNER_PARTICLE;
    }
  }
}


void setSourceTerm(PEM *pe, double dt){
  int i;
  double n0    = N0_forNumberDensity;
  double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    pe[i].sourceterm=0.0;
    if(pe[i].type==GHOST || pe[i].type== DUMMY_WALL ) continue;
    if(pe[i].boundary==INNER_PARTICLE){
      pe[i].sourceterm = gamma * (1.0/(dt*dt))*((pe[i].numdensity-n0)/n0);
    }else if(pe[i].boundary==SURFACE_PARTICLE){
      pe[i].sourceterm=0.0;
    }
  }
}


void setMatrix(PEM *pe, double dt){
  double xij, yij;
  double distance, distance2;
  double coefficientIJ;
  double n0 = N0_forLaplacian;
  int    i,j;
  double a;
  int n = NumberOfParticles;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    for(j=0;j<NumberOfParticles;j++){
      pe[i].coefficient[j] = 0.0;
    }
  }

  a = 2.0*DIM/(n0*Lambda);
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].boundary != INNER_PARTICLE) continue;
    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (pe[j].boundary==GHOST_OR_DUMMY) ) continue;
      xij = pe[j].x - pe[i].x;
      yij = pe[j].y - pe[i].y;
      distance2 = (xij*xij)+(yij*yij);
      distance  = sqrt(distance2);
      if(distance>=Re_forLaplacian)continue;
      coefficientIJ = a * weight(distance, Re_forLaplacian)/FluidDensity;
      pe[i].coefficient[j]  = (-1.0)*coefficientIJ;
      pe[i].coefficient[i] += coefficientIJ;
    }
    pe[i].coefficient[i] += (COMPRESSIBILITY)/(dt*dt);
  }
  exceptionalProcessingForBoundaryCondition(pe);
}


void exceptionalProcessingForBoundaryCondition(PEM *pe){
  /* If tere is no Dirichlet boundary condition on the fluid, 
     increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions. */
  checkBoundaryCondition(pe);
  increaseDiagonalTerm(pe);
}


void checkBoundaryCondition(PEM *pe){
  int i,j,count;
  double xij, yij, distance2;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if (pe[i].boundary==GHOST_OR_DUMMY){
      pe[i].boundaryflag=GHOST_OR_DUMMY;
    }else if (pe[i].boundary==SURFACE_PARTICLE){
      pe[i].boundaryflag=DIRICHLET_BOUNDARY_IS_CONNECTED;
    }else{
      pe[i].boundaryflag=DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
    }
  }

  do {
    count=0;
    for(i=0;i<NumberOfParticles;i++){
      if(pe[i].boundaryflag==DIRICHLET_BOUNDARY_IS_CONNECTED){
	for(j=0;j<NumberOfParticles;j++){
	  if( j==i ) continue;
	  if((pe[j].type==GHOST) || (pe[j].type== DUMMY_WALL)) continue;
	  if(pe[j].boundaryflag==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
	    xij = pe[j].x - pe[i].x;
	    yij = pe[j].y - pe[i].y;
	    distance2 = (xij*xij)+(yij*yij);
	    if(distance2>=Re2_forLaplacian)continue;
	    pe[j].boundaryflag=DIRICHLET_BOUNDARY_IS_CONNECTED;
	  }
	}
	pe[i].boundaryflag=DIRICHLET_BOUNDARY_IS_CHECKED;
	count++;
      }
    }
  } while (count!=0); /* This procedure is repeated until the all fluid or wall particles (which have Dirhchlet boundary condition in the particle group) are in the state of "DIRICHLET_BOUNDARY_IS_CHECKED".*/
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].boundaryflag==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
      fprintf(stderr,"WARNING: There is no dirichlet boundary condition for %d-th particle.\n",i );
    }
  }
}


void increaseDiagonalTerm(PEM *pe){
  int i;
  int n = NumberOfParticles;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<n;i++) {
    if(pe[i].boundaryflag == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED ){
      pe[i].coefficient[i] = 2.0 * pe[i].coefficient[i];
    }
  }
}


void solveSimultaniousEquationsByGaussEliminationMethod(PEM *pe){
  int    i,j,k;
  double c;
  double sumOfTerms;
  int    n = NumberOfParticles;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0; i<n; i++){ 
    pe[i].pressure = 0.0; 
  }
  #pragma omp parallel for  num_threads(8)
  for(i=0; i<n-1; i++){
    if ( pe[i].boundary != INNER_PARTICLE ) continue;
    for(j=i+1; j<n; j++){
      if(pe[j].boundary==GHOST_OR_DUMMY) continue;
      c = pe[j].coefficient[i]/pe[i].coefficient[i];
      for(k=i+1; k<n; k++){
	pe[j].coefficient[k] -= c * pe[i].coefficient[k];
      }
      pe[j].sourceterm -= c*pe[i].sourceterm;
    }
  }
  #pragma omp parallel for  num_threads(8)
  for( i=n-1; i>=0; i--){
    if ( pe[i].boundary != INNER_PARTICLE ) continue;
    sumOfTerms = 0.0;
    for( j=i+1; j<n; j++ ){
      if(pe[j].boundary==GHOST_OR_DUMMY) continue;
      sumOfTerms += pe[i].coefficient[j] * pe[j].pressure;
    }
    pe[i].pressure = (pe[i].sourceterm - sumOfTerms)/pe[i].coefficient[i];
  }
}


void removeNegativePressure(PEM *pe){
  int i;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++) {
    if(pe[i].pressure<0.0)pe[i].pressure=0.0;
  }
}


void setMinimumPressure(PEM *pe){
  double xij, yij, distance2;
  int i,j;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++) {
    if(pe[i].type==GHOST || pe[i].type==DUMMY_WALL)continue;
    pe[i].minpressure=pe[i].pressure;
    for(j=0;j<NumberOfParticles;j++) {
      if( (j==i) || (pe[j].type==GHOST) ) continue;
      if(pe[j].type==DUMMY_WALL) continue;
      xij = pe[j].x - pe[i].x;
      yij = pe[j].y - pe[i].y;
      distance2 = (xij*xij)+(yij*yij);
      if(distance2>=Re2_forGradient)continue;
      if( pe[i].minpressure > pe[j].pressure ){
	pe[i].minpressure = pe[j].pressure;
      }
    }
  }
}


void calPressureGradient(PEM *pe){
  int    i,j;
  double gradient_x, gradient_y;
  double xij, yij;
  double distance, distance2;
  double w,pij;
  double a;

  a =DIM/N0_forGradient;
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].type != FLUID) continue;
    gradient_x = 0.0;  gradient_y = 0.0;
    for(j=0;j<NumberOfParticles;j++){
      if( j==i ) continue;
      if( pe[j].type==GHOST ) continue;
      if( pe[j].type==DUMMY_WALL ) continue;
      xij = pe[j].x - pe[i].x;
      yij = pe[j].y - pe[i].y;
      distance2 = (xij*xij) + (yij*yij);
      distance = sqrt(distance2);
      if(distance<Re_forGradient){
	w =  weight(distance, Re_forGradient);
	pij = (pe[j].pressure - pe[i].minpressure)/distance2;
	gradient_x += xij*pij*w;
	gradient_y += yij*pij*w;
      }
    }
    gradient_x *= a;
    gradient_y *= a;
    pe[i].ax= (-1.0)*gradient_x/FluidDensity;
    pe[i].ay= (-1.0)*gradient_y/FluidDensity;
  }
}


void moveParticleUsingPressureGradient(PEM *pe, double dt){
  int i;
  
  #pragma omp parallel for  num_threads(8)
  for(i=0;i<NumberOfParticles;i++){
    if(pe[i].type == FLUID){
      pe[i].vx +=pe[i].ax*dt;
      pe[i].vy +=pe[i].ay*dt;

      pe[i].x +=pe[i].ax*dt*dt;
      pe[i].y +=pe[i].ay*dt*dt;
    }
    pe[i].ax=0.0;
    pe[i].ay=0.0;
  }
}

/*
void writeData_inProfFormat( void ){
  int i;
  FILE *fp;
  char fileName[256];

  sprintf(fileName, "output_%04d.prof",FileNumber);
  fp = fopen(fileName, "w");
  fprintf(fp,"%lf\n",Time);
  fprintf(fp,"%d\n",NumberOfParticles);
  for(i=0;i<NumberOfParticles;i++) {
    fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n"
	    ,pe[i].type, pe[i].x, pe[i].y, Position[i*3+2]
	    ,pe[i].vx, pe[i].vy, Velocity[i*3+2], pe[i].pressure, pe[i].numdensity);
  }
  fclose(fp);
  FileNumber++;
}


void writeData_inVtuFormat( void ){
  int i;
  double absoluteValueOfVelocity;
  FILE *fp;
  char fileName[1024];

  sprintf(fileName, "particle_%04d.vtu", FileNumber);
  fp=fopen(fileName,"w");
  fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",NumberOfParticles,NumberOfParticles);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%lf %lf %lf\n",pe[i].x,pe[i].y,Position[i*3+2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",pe[i].type);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    absoluteValueOfVelocity=
      sqrt( pe[i].vx*pe[i].vx + pe[i].vy*pe[i].vy + Velocity[i*3+2]*Velocity[i*3+2] );
    fprintf(fp,"%f\n",(float)absoluteValueOfVelocity);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)pe[i].pressure);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i+1);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='UInt8' Name='types' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"1\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
} 
*/

/*
//PEMの力の計算
void CalcForce(PEM *pe, double dt, double lc, double *param){
	int i,j;
	double cos_a, sin_a, ld; //衝突時の角度、中心間距離

	//2要素間の全ペアの力の計算
	for(i = 0; i < pe_n; i++){
		if(pe[i].exist == false) continue; //要素を使わない場合はforループを1周飛ばす

		//粒子-粒子間
		for(j = i+1; j < pe_n; j++){
			if(pe[j].exist == false) continue;
			if(pe[i].move == false && pe[j].move == false) continue; //両方動かない場合も力を計算する意味が無いので飛ばす。

			double lx = pe[j].x - pe[i].x;
			double ly = pe[j].y - pe[i].y;
			ld = sqrt(lx*lx+ly*ly);

			//ぶつかっている場合
			if(pe[i].r+pe[j].r > ld ){
				cos_a = lx/ld;
				sin_a = ly/ld;
				Force2Par(pe, pe_n, i, j, ld, cos_a, sin_a, dt);
			}else{
				pe[i].en[j] = 0.0;
				pe[i].es[j] = 0.0;
			}
		}
		//粒子-壁間だった場合
		for(j = pe_n; j < pe_n+wa_n; j++){
			if(pe[j].exist == false) continue;

			//直線の式
			//la*x + lb*y +lc = 0
			double la = sin(pe[j].phi);
			double lb = -cos(pe[j].phi);
			double lc = -(la*pe[j].x + lb*pe[j].y);
			double length_raw = la*pe[i].x + lb*pe[i].y + lc;
			ld = abs(length_raw);

			//ぶつかっている場合
			if(pe[i].r+pe[j].r - ld > 0 ){
				if(length_raw < 0){
					cos_a = la;
					sin_a = lb;
				}else{
					cos_a = -la;
					sin_a = -lb;
				}
				Force2Par(pe, pe_n, i, j, ld, cos_a, sin_a, dt);
			}else{
				pe[i].en[j] = 0.0;
				pe[i].es[j] = 0.0;
			}
		}
	}

	//外力の計算
	for(i = 0; i < pe_n+wa_n ; i++){
		if(pe[i].exist == false || pe[i].move == false) continue; //要素を使わないか、動かない場合は飛ばす
		pe[i].fy += -G * pe[i].m;   //重力
	}
}

//2要素の力を計算して代入する関数
void Force2Par(PEM* pe, int pe_n, int i, int j, double ld, double cos_a, double sin_a, double dt){
	//PEM用の係数変数
	//nが法線・sがせん断方向の値
	double un;	 double us;     //相対変位増分
	double veln; double vels;   //相対速度増分
	double hn;	 double hs;     //合力（法線,せん断）
	double kn;   double ks;
	double etan; double etas;
	double frc;

	//弾性・粘性係数の設定
	if(j<pe_n){
		kn   = 1000;    ks   = 5;
		etan = 50;     etas = 1;
		frc  = 10;
	}else{
		kn   = 50000;   ks   = 10000;
		etan = 10;     etas = 50;
		frc  = 1.0;
	}

	//法線方向成分とせん断方向成分の相対的変位増分の計算
	un   = +(pe[i].dx-pe[j].dx)*cos_a + (pe[i].dy-pe[j].dy)*sin_a;
	us   = -(pe[i].dx-pe[j].dx)*sin_a + (pe[i].dy-pe[j].dy)*cos_a 
		+(pe[i].r*pe[i].dphi + pe[j].r*pe[j].dphi);

	veln = +(pe[i].vx-pe[j].vx)*cos_a + (pe[i].vy-pe[j].vy)*sin_a;
	vels = -(pe[i].vx-pe[j].vx)*sin_a + (pe[i].vy-pe[j].vy)*cos_a
		+(pe[i].r*pe[i].vphi + pe[j].r*pe[j].vphi);

	//法線方向成分とせん断方向成分の合力の計算
	//バネの力を追加
	pe[i].en[j] = pe[i].en[j] + kn * un;
	pe[i].es[j] = pe[i].es[j] + ks * us;
	//ダンパの力を追加
	hn = pe[i].en[j] + etan * veln;
	hs = pe[i].es[j] + etas * vels;

	if(hn <= 0.0){
		hs=0.0;
	}else if((abs(hs)-frc*hn) >= 0.0){
		hs=frc*fabs(hn)*hs/fabs(hs);
	}

	//粒子i,粒子jのx方向成分とy方向成分の合力の計算
	pe[i].fx = -hn*cos_a + hs*sin_a + pe[i].fx;
	pe[i].fy = -hn*sin_a - hs*cos_a + pe[i].fy;
	pe[i].fm = pe[i].fm - pe[i].r*hs;

	pe[j].fx = hn*cos_a - hs*sin_a + pe[j].fx;
	pe[j].fy = hn*sin_a + hs*cos_a + pe[j].fy;
	pe[j].fm = pe[j].fm - pe[j].r*hs;
}
*/
//おまじない
#pragma managed
