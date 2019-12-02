//���܂��Ȃ�
//(�Ӗ�: ���̃t�@�C��������#include���Ă�����������ǂ܂Ȃ�����)
#ifndef __PEMSTRUCT_H__
#define __PEMSTRUCT_H__


//�v�f1�̏��������Ă����\����
//���ɑ������~�����Ȃ����炱���ɒǉ�
struct PEM{
	
	bool exist;				//�v�f���g���Ă��邩�ǂ���
	bool move;				//�v�f���������ǂ���
	int type;				//���q�̃^�C�v�ԍ�

	double r;				//���q�̔��a
	double m;				//����
	double Ir;				//�������[�����g
	double x, y, phi;		//���W(x,y),�p�xphi
	double vx, vy, vphi;	//���x
	double ax, ay, aphi;	//�����x
	double dx, dy, dphi;    //�ψʑ����Ɖ�]�ψʑ���
	double fx, fy, fm;		//���͂ƃ��[�����g
	double *en, *es;        //���q�̖@�������e���͂Ƃ���f�����e����
	double numdensity;		//���q�����x
	int boundary;			//�f�B���N�����E������t�����邩�ǂ����̃t���O
	double sourceterm;		//���͂̃|�A�\���������̉E�Ӄx�N�g��
	double *coefficient;		//�A���ꎟ�������̌W���s��
	int boundaryflag;		//���E�������t������Ă��邩�ǂ����̃t���O
	double pressure;		//����
	double minpressure;		//���q�̎���̍ŏ�����
};

//Element.type���������Ƃ킩��ɂ����̂ŕʖ������Ă���

//���܂��Ȃ�
#endif
