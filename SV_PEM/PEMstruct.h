//おまじない
//(意味: このファイルが既に#includeしてあったらもう読まないこと)
#ifndef __PEMSTRUCT_H__
#define __PEMSTRUCT_H__


//要素1つの情報を書いておく構造体
//他に属性が欲しくなったらここに追加
struct PEM{
	
	bool exist;				//要素が使われているかどうか
	bool move;				//要素が動くかどうか
	int type;				//粒子のタイプ番号

	double r;				//粒子の半径
	double m;				//質量
	double Ir;				//慣性モーメント
	double x, y, phi;		//座標(x,y),角度phi
	double vx, vy, vphi;	//速度
	double ax, ay, aphi;	//加速度
	double dx, dy, dphi;    //変位増分と回転変位増分
	double fx, fy, fm;		//合力とモーメント
	double *en, *es;        //粒子の法線方向弾性力とせん断方向弾性力
	double numdensity;		//粒子数密度
	int boundary;			//ディリクレ境界条件を付加するかどうかのフラグ
	double sourceterm;		//圧力のポアソン方程式の右辺ベクトル
	double *coefficient;		//連立一次方程式の係数行列
	int boundaryflag;		//境界条件が付加されているかどうかのフラグ
	double pressure;		//圧力
	double minpressure;		//粒子の周りの最小圧力
};

//Element.typeが数字だとわかりにくいので別名をつけておく

//おまじない
#endif
