//おまじない
//(意味: ここから下は普通のC++の構文で書いてありますよ)
#pragma unmanaged

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PEMstruct.h"

const double PI = 3.1415926535;		//円周率
const double G = 9.80665;		    //重力加速度
const double rho = 10;		        //粒子の密度

//力の計算をする関数
void CalcForce(int pe_n, int wa_n, PEM *pe, double dt, double lc);

//PEMに値を一気にセットする便利関数
//PEMを書き換えたらここも変えておくとよい
//MD用セット関数(MDに用いない回転方向のパラメータは0.0)
void setParam(PEM *pe, int i, int n, bool move, double r, double x, double y, double vx, double vy, double m){
	pe[i].exist = true; //計算で使うよう設定

	//引数で与えられた値の設定
	pe[i].move = move;
	pe[i].r = r;//表示上の大きさ
	pe[i].x = x;   pe[i].y = y;   pe[i].phi = 0.0;
	pe[i].vx = vx; pe[i].vy = vy; pe[i].vphi = 0.0;

	//加速度,変位増分,力を0で初期化
	pe[i].ax = pe[i].ay = pe[i].aphi = 0.0;
	pe[i].dx = pe[i].dy = pe[i].dphi = 0.0;
	pe[i].fx = pe[i].fy = pe[i].fm   = 0.0;

	//質量,慣性モーメントの設定
	pe[i].m = m;
	pe[i].Ir = 0.0;

	//弾性力の初期化
	for(int j=0; j<n ;j++){
		pe[i].en[j] = 0.0;
		pe[i].es[j] = 0.0;
	}
}

//系の中の各要素の情報を細かく設定する関数
//この関数を呼ぶと粒子とかが初期状態になる
void Config(int pe_n, int wa_n, PEM *pe){

	int n=pe_n+wa_n;

	//存在初期化
	for(int i=0 ;i<n ;i++){
		pe[i].exist=false;
	}

	//粒子の初期条件
	//setParamMD(pe,  i, n, move,          r,        x,        y, vx, vy,             m);
	setParam(  pe,  0, n, true,  3.0*1e-11,  1*1e-10, -1*1e-10,  0,  0, 27.0*1.67e-27);
	setParam(  pe,  1, n, true,  3.0*1e-11,  1*1e-10,  2*1e-10,  0,  0, 27.0*1.67e-27);
	setParam(  pe,  2, n, true,  3.0*1e-11, -1*1e-10,  0*1e-10,  0,  0, 27.0*1.67e-27);
}

//1ステップの動きを計算する関数
void CalcStep(int pe_n, int wa_n, PEM *pe, double dt, double lc, double *param){
	int i;

	//各粒子の力の初期化
	for(i = 0; i < pe_n+wa_n; i++){
		pe[i].fx=0.0; pe[i].fy=0.0; pe[i].fm=0.0;
	}

	//力の計算
	CalcForce(pe_n, wa_n, pe, dt, lc);

	//-----------------ここから位置などの計算------------

	//力から加速度，速度，変位の更新
	for(i = 0; i < pe_n+wa_n; i++){
		if(pe[i].exist == false || pe[i].move == false) continue; //要素を使わないか、動かない場合は飛ばす

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

		//周期境界条件による位置の更新
		pe[i].x=pe[i].x-floor((pe[i].x+lc/2.0)/lc)*lc;
		pe[i].y=pe[i].y-floor((pe[i].y+lc/2.0)/lc)*lc;

		//phiの範囲を[0, 2PI]に収める(収めなくても問題なく動く)
		if(pe[i].phi < 0){
			pe[i].phi += 2*PI;
		}else if(2*PI < pe[i].phi){
			pe[i].phi -= 2*PI;
		}

		//エネルギーなどの計算結果
		param[1]=0.0;
	}
}

//MDの力の計算
void CalcForce(int pe_n, int wa_n, PEM *pe, double dt, double lc){
	int i,j;

	//モースポテンシャル関連定数
	double dkb    = 1.38062e-23;  //ボルツマン定数
	double ev     = 1.60219e-19;  //電荷
	double epsilon= 0.27*ev;	  //イプシロン
	double alpha  = 1.16e+10;     //α
	double r0     = 3.25e-10;     //r0

	//モースポテンシャル原子間力計算用変数
	double sbr,sbr2;
	double fd,fda,fdb;	
	double fdx1,fdy1;	
	double fdx2,fdy2;
	double rx_ij,rx_ji;	//原子間変位ij，ji(x方向)
	double ry_ij,ry_ji;	//原子間変位ij，ji(y方向)
	double rd;

	//2要素間の全ペアの力の計算
	for(i = 0; i < pe_n; i++){
		if(pe[i].exist == false) continue; //要素を使わない場合はforループを1周飛ばす

		//粒子-粒子間
		for(j = i+1; j < pe_n; j++){
			if(pe[j].exist == false) continue;
			if(pe[i].move == false && pe[j].move == false) continue; //両方動かない場合も力を計算する意味が無いので飛ばす。

			//周期境界条件による粒子間距離の計算とカットオフ
			rx_ij=pe[i].x-pe[j].x;
			ry_ij=pe[i].y-pe[j].y;
			rx_ij=rx_ij-floor(rx_ij/lc+0.5)*lc;
			ry_ij=ry_ij-floor(ry_ij/lc+0.5)*lc;
			rx_ji=-rx_ij;
			ry_ji=-ry_ij;
			rd = sqrt(rx_ij*rx_ij+ry_ij*ry_ij);

			//原子間の相互作用の計算
			if(rd>=lc/2.0){
				fd=0.0;
			}else{
				sbr = r0-rd;
				sbr2 = alpha*sbr;
				fda = exp(sbr2);
				fdb = exp(2.0*sbr2);
				fd=(-2.0*fdb+fda)*alpha*epsilon;
			}

			//作用力の計算
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

//おまじない
#pragma managed
