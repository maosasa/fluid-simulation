#pragma once

#include <stdio.h>
#include <math.h>
#include "PEMstruct.h"	//���q�p�����[�^�����\���̂��L�q����Ă���

#define ITERATE_NUM 10
#define LC 1.0

void Config(int rows, int pe_in_raw, int wa_n, PEM *pe, double kv);
void CalcStep(int pe_n, int wa_n, PEM *pe, double dt, double lc, double *param);

namespace SV_PEM {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace System::Drawing::Drawing2D;

	public ref class Form1 : public System::Windows::Forms::Form
	{
	private:int rows;
	private:int kv;				//���S���W��
	private:int pe_in_row;
	private:int pe_n;			//���q��
	private:int wa_n;			//�ǂ̐�
	private:double lc;          //��\�X�P�[��(�Z���̑傫��(MD)/�ǂ̈ʒu(DEM))
	private:int step;			//�v�Z�X�e�b�v��
	private:double dt;			//�v�Z�̎��ԍ���(s)
	private:PEM *pe;            //���q�p�����[�^�����\���̃|�C���^
	private:double camscale;	//�J�����{��
	private:double *param;		//�v�Z���ʂ̎󂯓n���p�ϐ��|�C���^

	private:Graphics^ gr;		//�`��p�|�C���^1
	private:Graphics^ gr2;		//�`��p�|�C���^2
	private:Pen^ blackPen;		//�`���p�|�C���^�i���j
	private:SolidBrush^ redBrush;	    //�h�ׂ��p�|�C���^�i�ԁj
	private:SolidBrush^ blackBrush;	    //�h�ׂ��p�|�C���^�i���j
	private:SolidBrush^ blueBrush;	    //�h�ׂ��p�|�C���^�i�j
	private:SolidBrush^ greenBrush;	    //�h�ׂ��p�|�C���^�i�΁j
	private:array<SolidBrush^>^ atomBrush;	        //�e���q�̐F
	private:BufferedGraphicsContext^ myContext;		//�摜�o�b�t�@
	private:BufferedGraphics^ myBuffer;				//�摜�o�b�t�@
	
	private: System::Windows::Forms::Label^  s;
	private: System::Windows::Forms::Label^  scale;
	private: System::Windows::Forms::Label^  time;
	private: System::Windows::Forms::Label^  ScaleLabel;
	private: System::Windows::Forms::Button^  BiggerButton;
	private: System::Windows::Forms::Button^  SmallerButton;
	private: System::Windows::Forms::Label^  TimeLabel;
	private: System::Windows::Forms::Button^  ResetButton;
	private: System::Windows::Forms::TrackBar^  trackBar1;
	private: System::Windows::Forms::Label^  label1;

	private: System::Windows::Forms::PictureBox^  pictureBox1;	

	public:
		/// <summary>
		/// �R���X�g���N�^
		/// �N���X�������ɍŏ��ɌĂ΂��
		/// </summary>
		Form1()
		{
			InitializeComponent();
			//dt=1e-15;                 //�v�Z�̎��ԍ���(s)(MD)
			dt=0.005;					//�v�Z�̎��ԍ���(s)(DEM)
			//lc=1e-9;				    //MD�̑�\�T�C�Y
			lc=LC;					    //DEM�̑�\�T�C�Y
			kv=0.1;
			pe_n=5000;					//���q���̐ݒ�
			wa_n=4;						//�ǂ̐��̐ݒ�
			param = new double[10];		//�v�Z���ʂȂǂ̊i�[�p�̃������m��
			pe = new PEM[pe_n+wa_n];	//���q�\���̂̃������m��
			for(int i=0; i<pe_n+wa_n ;i++){
				pe[i].en = new double[pe_n+wa_n];
				pe[i].es = new double[pe_n+wa_n];
				pe[i].coefficient = new double[pe_n+wa_n];
			}
			Config(rows, pe_in_row, wa_n, pe, kv); //���q�̏����l�ݒ�(main_MD.cpp�Q��)

			atomBrush = gcnew array<SolidBrush^>(pe_n);     //���q�̐F
			redBrush = gcnew SolidBrush( Color::Red );		//�h�ׂ��i�ԁj
			blackBrush = gcnew SolidBrush( Color::Black );	//�h�ׂ��i���j
			blueBrush = gcnew SolidBrush( Color::Blue );	//�h�ׂ��i�j
			greenBrush = gcnew SolidBrush( Color::Green );	//�h�ׂ��i�j
			blackPen = gcnew Pen(System::Drawing::Color::Black, 3);	//�`���i���j
			gr = pictureBox1->CreateGraphics();						//�`��
		}

	protected:
		/// <summary>
		/// �f�X�g���N�^
		/// �g�p���̃��\�[�X�����ׂăN���[���A�b�v����
		/// </summary>
		~Form1()
		{
			if (components)
			{
				delete components;
			}
			delete[] pe;
		}

	private: System::Windows::Forms::Timer^  timer1;
	private: System::Windows::Forms::Button^  StopButton;
	private: System::Windows::Forms::Button^  StartButton;
	private: System::ComponentModel::IContainer^  components;
	private:

#pragma region Windows Form Designer generated code
		/// <summary>
		/// �f�U�C�i �T�|�[�g�ɕK�v�ȃ��\�b�h�ł��B���̃��\�b�h�̓��e��
		/// �R�[�h �G�f�B�^�ŕύX���Ȃ��ł��������B
		/// ����� - ���������ƂŃR���p�N�g�ɂł��܂��B
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			this->timer1 = (gcnew System::Windows::Forms::Timer(this->components));
			this->StopButton = (gcnew System::Windows::Forms::Button());
			this->StartButton = (gcnew System::Windows::Forms::Button());
			this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			this->s = (gcnew System::Windows::Forms::Label());
			this->scale = (gcnew System::Windows::Forms::Label());
			this->time = (gcnew System::Windows::Forms::Label());
			this->ScaleLabel = (gcnew System::Windows::Forms::Label());
			this->BiggerButton = (gcnew System::Windows::Forms::Button());
			this->SmallerButton = (gcnew System::Windows::Forms::Button());
			this->TimeLabel = (gcnew System::Windows::Forms::Label());
			this->ResetButton = (gcnew System::Windows::Forms::Button());
			this->trackBar1 = (gcnew System::Windows::Forms::TrackBar());
			this->label1 = (gcnew System::Windows::Forms::Label());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pictureBox1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->trackBar1))->BeginInit();
			this->SuspendLayout();
			// 
			// timer1
			// 
			this->timer1->Tick += gcnew System::EventHandler(this, &Form1::timer1_Tick);
			// 
			// StopButton
			// 
			this->StopButton->Location = System::Drawing::Point(400, 483);
			this->StopButton->Name = L"StopButton";
			this->StopButton->Size = System::Drawing::Size(75, 23);
			this->StopButton->TabIndex = 2;
			this->StopButton->Text = L"STOP";
			this->StopButton->UseVisualStyleBackColor = true;
			this->StopButton->Click += gcnew System::EventHandler(this, &Form1::StopButton_Click);
			// 
			// StartButton
			// 
			this->StartButton->Location = System::Drawing::Point(308, 483);
			this->StartButton->Name = L"StartButton";
			this->StartButton->Size = System::Drawing::Size(75, 23);
			this->StartButton->TabIndex = 3;
			this->StartButton->Text = L"START";
			this->StartButton->UseVisualStyleBackColor = true;
			this->StartButton->Click += gcnew System::EventHandler(this, &Form1::StartButton_Click);
			// 
			// pictureBox1
			// 
			this->pictureBox1->BackColor = System::Drawing::SystemColors::Window;
			this->pictureBox1->Location = System::Drawing::Point(25, 21);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(450, 450);
			this->pictureBox1->TabIndex = 4;
			this->pictureBox1->TabStop = false;
			// 
			// s
			// 
			this->s->AutoSize = true;
			this->s->Font = (gcnew System::Drawing::Font(L"MS UI Gothic", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(128)));
			this->s->Location = System::Drawing::Point(99, 512);
			this->s->Name = L"s";
			this->s->Size = System::Drawing::Size(25, 16);
			this->s->TabIndex = 7;
			this->s->Text = L"[s]";
			// 
			// scale
			// 
			this->scale->AutoSize = true;
			this->scale->Location = System::Drawing::Point(132, 492);
			this->scale->Name = L"scale";
			this->scale->Size = System::Drawing::Size(32, 12);
			this->scale->TabIndex = 8;
			this->scale->Text = L"scale";
			// 
			// time
			// 
			this->time->AutoSize = true;
			this->time->Location = System::Drawing::Point(29, 492);
			this->time->Name = L"time";
			this->time->Size = System::Drawing::Size(27, 12);
			this->time->TabIndex = 9;
			this->time->Text = L"time";
			// 
			// ScaleLabel
			// 
			this->ScaleLabel->AutoSize = true;
			this->ScaleLabel->Font = (gcnew System::Drawing::Font(L"MS UI Gothic", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(128)));
			this->ScaleLabel->Location = System::Drawing::Point(130, 512);
			this->ScaleLabel->Name = L"ScaleLabel";
			this->ScaleLabel->Size = System::Drawing::Size(16, 16);
			this->ScaleLabel->TabIndex = 10;
			this->ScaleLabel->Text = L"1";
			// 
			// BiggerButton
			// 
			this->BiggerButton->Location = System::Drawing::Point(418, 512);
			this->BiggerButton->Name = L"BiggerButton";
			this->BiggerButton->Size = System::Drawing::Size(20, 23);
			this->BiggerButton->TabIndex = 11;
			this->BiggerButton->Text = L"<";
			this->BiggerButton->UseVisualStyleBackColor = true;
			this->BiggerButton->Click += gcnew System::EventHandler(this, &Form1::BiggerButton_Click);
			// 
			// SmallerButton
			// 
			this->SmallerButton->Location = System::Drawing::Point(444, 512);
			this->SmallerButton->Name = L"SmallerButton";
			this->SmallerButton->Size = System::Drawing::Size(20, 23);
			this->SmallerButton->TabIndex = 12;
			this->SmallerButton->Text = L">";
			this->SmallerButton->UseVisualStyleBackColor = true;
			this->SmallerButton->Click += gcnew System::EventHandler(this, &Form1::SmallerButton_Click);
			// 
			// TimeLabel
			// 
			this->TimeLabel->AutoSize = true;
			this->TimeLabel->Font = (gcnew System::Drawing::Font(L"MS UI Gothic", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(128)));
			this->TimeLabel->Location = System::Drawing::Point(28, 512);
			this->TimeLabel->Name = L"TimeLabel";
			this->TimeLabel->Size = System::Drawing::Size(16, 16);
			this->TimeLabel->TabIndex = 13;
			this->TimeLabel->Text = L"0";
			// 
			// ResetButton
			// 
			this->ResetButton->Location = System::Drawing::Point(308, 512);
			this->ResetButton->Name = L"ResetButton";
			this->ResetButton->Size = System::Drawing::Size(75, 23);
			this->ResetButton->TabIndex = 14;
			this->ResetButton->Text = L"RESET";
			this->ResetButton->UseVisualStyleBackColor = true;
			this->ResetButton->Click += gcnew System::EventHandler(this, &Form1::ResetButton_Click);
			// 
			// trackBar1
			// 
			this->trackBar1->Location = System::Drawing::Point(198, 492);
			this->trackBar1->Name = L"trackBar1";
			this->trackBar1->Size = System::Drawing::Size(104, 45);
			this->trackBar1->TabIndex = 15;
			this->trackBar1->Scroll += gcnew System::EventHandler(this, &Form1::trackBar1_Scroll);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(228, 540);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(35, 12);
			this->label1->TabIndex = 16;
			this->label1->Text = L"label1";
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(504, 572);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->trackBar1);
			this->Controls->Add(this->ResetButton);
			this->Controls->Add(this->TimeLabel);
			this->Controls->Add(this->SmallerButton);
			this->Controls->Add(this->BiggerButton);
			this->Controls->Add(this->ScaleLabel);
			this->Controls->Add(this->time);
			this->Controls->Add(this->scale);
			this->Controls->Add(this->s);
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->StartButton);
			this->Controls->Add(this->StopButton);
			this->Name = L"Form1";
			this->Text = L"Form1";
			this->Load += gcnew System::EventHandler(this, &Form1::Form1_Load);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->pictureBox1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^  >(this->trackBar1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

			//�v���O�������N�����ăE�B���h�E���`�悳��钼�O�ɌĂ΂��֐�
			//�����ݒ�Ȃǂ��Ă����Ƃ悢
	private: System::Void Form1_Load(System::Object^  sender, System::EventArgs^  e) {
				 //�����l
				 timer1->Interval = 1;     // ���ؕb���ɲ���Ĕ��������邩�w��
				 timer1->Enabled = false;  // ��ϰ���~��Ԃŏ�����
				 step=0;
				 //�J�����X�P�[��(=Form�̕�����\�T�C�Y)
				 camscale=(double)pictureBox1->Width/lc;
				 ScaleLabel->Text = Convert::ToString(camscale);
			 }

			 //START�{�^�����N���b�N���ꂽ�Ƃ��ȉ������s
	private: System::Void StartButton_Click(System::Object^  sender, System::EventArgs^  e) {	
				 label1->Text=trackBar1->Value.ToString("F8");
				 kv = 10^(-1*(trackBar1->Value));
				 Config(rows, pe_in_row, wa_n, pe, kv);
				 timer1->Start();


			 }

			 //STOP�{�^�����N���b�N���ꂽ�Ƃ��ȉ������s
	private: System::Void StopButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 timer1->Stop();
			 }

			 //�X�P�[���g��
	private: System::Void BiggerButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 camscale *= 1.5;
				 ScaleLabel->Text = Convert::ToString(camscale);
			 }

			 //�X�P�[���k��
	private: System::Void SmallerButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 camscale /= 1.5;
				 ScaleLabel->Text = Convert::ToString(camscale);
			 }

			 //���������̍ēǂݍ���(RESET)
	private: System::Void ResetButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 for(int i=0; i<pe_n+wa_n ;i++){
					 pe[i].en = new double[pe_n+wa_n];
					 pe[i].es = new double[pe_n+wa_n];
				 }
				 Config(rows, pe_in_row, wa_n, pe, kv); //���q�̏����l�ݒ�(PEM.cpp�Q��)

				 timer1->Interval = 1;     // ���ؕb���ɲ���Ĕ��������邩�w��
				 timer1->Enabled = false;  // ��ϰ���~��Ԃŏ�����
				 step=0;
				 //�J�����X�P�[��(=Form�̕�����\�T�C�Y)
				 camscale=(double)pictureBox1->Width/lc;
				 ScaleLabel->Text = Convert::ToString(camscale);
			 }

			 //�^�C�}�[��Tick�Ōv�Z����
	private: System::Void timer1_Tick(System::Object^  sender, System::EventArgs^  e) {

				 //���q�@(PEM)�̃��C���v���O����
				 //1Tick��10��v�Z����
				 for(int i = 0; i < ITERATE_NUM; i++){
					 CalcStep(pe_n, wa_n, pe, dt, lc, param);
					 step++;
				 }

				 //param�ɏ������܂ꂽ�G�l���M�[�Ȃǂ̏���Form1�ɕ\��


				 /// <summary>
				 /// �ȉ��A��ʂւ̕`��R�[�h
				 /// </summary>

				 myContext = gcnew BufferedGraphicsContext();	
				 myBuffer = myContext->Allocate(pictureBox1->CreateGraphics(),pictureBox1->DisplayRectangle);
				 myBuffer->Graphics->Clear(Color::White);
				 double camx = (double)pictureBox1->Width;	//pictureBox�̕�
				 double camy = (double)pictureBox1->Height;	//pictureBox�̍���

				 //���q�̏ꍇ
				 for(int ni=0; ni < pe_n+wa_n; ni++){
						 // �\�����q�T�C�Y(��)
						 float radius = (float)(camscale*pe[ni].r);
						 // �`�悷�闱�q�̈ʒu�̎w��
						 float xd = (float)(camx/2+camscale*pe[ni].x);
						 float yd = (float)(camy/2-camscale*pe[ni].y);

					 if(pe[ni].type==0){
						 //myBuffer->Graphics->DrawEllipse(blackPen,xd-radius,yd-radius,2*radius,2*radius);
						 myBuffer->Graphics->FillEllipse(blueBrush,xd-radius,yd-radius,2*radius,2*radius);
					 }else{
						 
						 myBuffer->Graphics->FillEllipse(blackBrush,xd-radius,yd-radius,2*radius,2*radius);
					 }
				 }
				 
				 //�ǂ̏ꍇ
				 for(int ni=pe_n; ni < pe_n+wa_n; ni++){
					 //�ǂ��Ȃ����Ƃɐݒ肳��Ă������΂�
					 if(pe[ni].exist == false) continue;

					 //�����̎�
					 //la*x + lb*y + lc = 0
					 double la =  sin(pe[ni].phi);
					 double lb = -cos(pe[ni].phi);
					 double lc = -(la*pe[ni].x + lb*pe[ni].y);
					 /*
					 //���̗��[�̃V�~�����[�V������̍��W
					 double x1, y1, x2, y2;

					 if(fabs(lb) < 0.01){ //y���ɂ����������s�ȏꍇ
						 y1 = -camx/2.0 / camscale;
						 y2 = camx/2.0 / camscale;

						 x1 = (-lc-lb*y1)/la;
						 x2 = (-lc-lb*y2)/la;
					 }else{
						 x1 = -camx/2.0 / camscale; //x�������[�_�̌����W�ł̒l
						 x2 = camx/2.0 / camscale;

						 y1 = (-lc-la*x1)/lb;
						 y2 = (-lc-la*x2)/lb;
					 }
					 
					 //�J�������W�ɂ���
					 float camx1 = (float)(camx/2.0 + camscale*x1);
					 float camy1 = (float)(camx/2.0 + camscale*y1);
					 float camx2 = (float)(camx/2.0 + camscale*x2);
					 float camy2 = (float)(camx/2.0 + camscale*y2);

					 //��
					 float camwidth = (float)(camscale * 2.0 * pe[ni].r);

					 //�y���ɕ���ݒ�
					 //C++�̕��@�ɂȂ��������B���܂�C�ɂ��Ȃ��Ă悢
					 Pen^ linePen = gcnew Pen(System::Drawing::Color::Black, camwidth);
					 //���̕`��
					 myBuffer->Graphics->DrawLine(linePen, camx1, camy1, camx2, camy2);
					 */
				 }
				 
				 // ���Ԃ̕\��
				 double tt = step*dt;
				 TimeLabel->Text=Convert::ToString(tt);
				 myBuffer->Render(gr);
			 }
	private: System::Void trackBar1_Scroll(System::Object^  sender, System::EventArgs^  e) {
				 
			 }
};
}

