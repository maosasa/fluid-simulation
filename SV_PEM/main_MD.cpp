//���܂��Ȃ�
//(�Ӗ�: �������牺�͕��ʂ�C++�̍\���ŏ����Ă���܂���)
#pragma unmanaged

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PEMstruct.h"

const double PI = 3.1415926535;		//�~����
const double G = 9.80665;		    //�d�͉����x
const double rho = 10;		        //���q�̖��x

//�͂̌v�Z������֐�
void CalcForce(int pe_n, int wa_n, PEM *pe, double dt, double lc);

//PEM�ɒl����C�ɃZ�b�g����֗��֐�
//PEM�������������炱�����ς��Ă����Ƃ悢
//MD�p�Z�b�g�֐�(MD�ɗp���Ȃ���]�����̃p�����[�^��0.0)
void setParam(PEM *pe, int i, int n, bool move, double r, double x, double y, double vx, double vy, double m){
	pe[i].exist = true; //�v�Z�Ŏg���悤�ݒ�

	//�����ŗ^����ꂽ�l�̐ݒ�
	pe[i].move = move;
	pe[i].r = r;//�\����̑傫��
	pe[i].x = x;   pe[i].y = y;   pe[i].phi = 0.0;
	pe[i].vx = vx; pe[i].vy = vy; pe[i].vphi = 0.0;

	//�����x,�ψʑ���,�͂�0�ŏ�����
	pe[i].ax = pe[i].ay = pe[i].aphi = 0.0;
	pe[i].dx = pe[i].dy = pe[i].dphi = 0.0;
	pe[i].fx = pe[i].fy = pe[i].fm   = 0.0;

	//����,�������[�����g�̐ݒ�
	pe[i].m = m;
	pe[i].Ir = 0.0;

	//�e���͂̏�����
	for(int j=0; j<n ;j++){
		pe[i].en[j] = 0.0;
		pe[i].es[j] = 0.0;
	}
}

//�n�̒��̊e�v�f�̏����ׂ����ݒ肷��֐�
//���̊֐����ĂԂƗ��q�Ƃ���������ԂɂȂ�
void Config(int pe_n, int wa_n, PEM *pe){

	int n=pe_n+wa_n;

	//���ݏ�����
	for(int i=0 ;i<n ;i++){
		pe[i].exist=false;
	}

	//���q�̏�������
	//setParamMD(pe,  i, n, move,          r,        x,        y, vx, vy,             m);
	setParam(  pe,  0, n, true,  3.0*1e-11,  1*1e-10, -1*1e-10,  0,  0, 27.0*1.67e-27);
	setParam(  pe,  1, n, true,  3.0*1e-11,  1*1e-10,  2*1e-10,  0,  0, 27.0*1.67e-27);
	setParam(  pe,  2, n, true,  3.0*1e-11, -1*1e-10,  0*1e-10,  0,  0, 27.0*1.67e-27);
}

//1�X�e�b�v�̓������v�Z����֐�
void CalcStep(int pe_n, int wa_n, PEM *pe, double dt, double lc, double *param){
	int i;

	//�e���q�̗͂̏�����
	for(i = 0; i < pe_n+wa_n; i++){
		pe[i].fx=0.0; pe[i].fy=0.0; pe[i].fm=0.0;
	}

	//�͂̌v�Z
	CalcForce(pe_n, wa_n, pe, dt, lc);

	//-----------------��������ʒu�Ȃǂ̌v�Z------------

	//�͂�������x�C���x�C�ψʂ̍X�V
	for(i = 0; i < pe_n+wa_n; i++){
		if(pe[i].exist == false || pe[i].move == false) continue; //�v�f���g��Ȃ����A�����Ȃ��ꍇ�͔�΂�

		pe[i].ax   = pe[i].fx / pe[i].m;
		pe[i].ay   = pe[i].fy / pe[i].m;
		pe[i].aphi = pe[i].fm / pe[i].Ir;

		pe[i].vx   += dt * pe[i].ax;
		pe[i].vy   += dt * pe[i].ay;
		pe[i].vphi += dt * pe[i].aphi;

		pe[i].dx   = dt * pe[i].vx;
		pe[i].dy   = dt * pe[i].vy;
		pe[i].dphi = dt * pe[i].vphi;

		pe[i].x   += pe[i].dx;
		pe[i].y   += pe[i].dy;
		pe[i].phi += pe[i].dphi;

		//�������E�����ɂ��ʒu�̍X�V
		pe[i].x=pe[i].x-floor((pe[i].x+lc/2.0)/lc)*lc;
		pe[i].y=pe[i].y-floor((pe[i].y+lc/2.0)/lc)*lc;

		//phi�͈̔͂�[0, 2PI]�Ɏ��߂�(���߂Ȃ��Ă����Ȃ�����)
		if(pe[i].phi < 0){
			pe[i].phi += 2*PI;
		}else if(2*PI < pe[i].phi){
			pe[i].phi -= 2*PI;
		}

		//�G�l���M�[�Ȃǂ̌v�Z����
		param[1]=0.0;
	}
}

//MD�̗͂̌v�Z
void CalcForce(int pe_n, int wa_n, PEM *pe, double dt, double lc){
	int i,j;

	//���[�X�|�e���V�����֘A�萔
	double dkb    = 1.38062e-23;  //�{���c�}���萔
	double ev     = 1.60219e-19;  //�d��
	double epsilon= 0.27*ev;	  //�C�v�V����
	double alpha  = 1.16e+10;     //��
	double r0     = 3.25e-10;     //r0

	//���[�X�|�e���V�������q�ԗ͌v�Z�p�ϐ�
	double sbr,sbr2;
	double fd,fda,fdb;	
	double fdx1,fdy1;	
	double fdx2,fdy2;
	double rx_ij,rx_ji;	//���q�ԕψ�ij�Cji(x����)
	double ry_ij,ry_ji;	//���q�ԕψ�ij�Cji(y����)
	double rd;

	//2�v�f�Ԃ̑S�y�A�̗͂̌v�Z
	for(i = 0; i < pe_n; i++){
		if(pe[i].exist == false) continue; //�v�f���g��Ȃ��ꍇ��for���[�v��1����΂�

		//���q-���q��
		for(j = i+1; j < pe_n; j++){
			if(pe[j].exist == false) continue;
			if(pe[i].move == false && pe[j].move == false) continue; //���������Ȃ��ꍇ���͂��v�Z����Ӗ��������̂Ŕ�΂��B

			//�������E�����ɂ�闱�q�ԋ����̌v�Z�ƃJ�b�g�I�t
			rx_ij=pe[i].x-pe[j].x;
			ry_ij=pe[i].y-pe[j].y;
			rx_ij=rx_ij-floor(rx_ij/lc+0.5)*lc;
			ry_ij=ry_ij-floor(ry_ij/lc+0.5)*lc;
			rx_ji=-rx_ij;
			ry_ji=-ry_ij;
			rd = sqrt(rx_ij*rx_ij+ry_ij*ry_ij);

			//���q�Ԃ̑��ݍ�p�̌v�Z
			if(rd>=lc/2.0){
				fd=0.0;
			}else{
				sbr = r0-rd;
				sbr2 = alpha*sbr;
				fda = exp(sbr2);
				fdb = exp(2.0*sbr2);
				fd=(-2.0*fdb+fda)*alpha*epsilon;
			}

			//��p�͂̌v�Z
			fdx1=-fd*rx_ij/rd;
			fdy1=-fd*ry_ij/rd;
			fdx2=-fd*rx_ji/rd;
			fdy2=-fd*ry_ji/rd;

			pe[i].fx=pe[i].fx+fdx1;
			pe[i].fy=pe[i].fy+fdy1;
			pe[j].fx=pe[j].fx+fdx2;
			pe[j].fy=pe[j].fy+fdy2;
		}
	}
}

//���܂��Ȃ�
#pragma managed
