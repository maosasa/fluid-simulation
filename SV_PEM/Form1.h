#pragma once

#include <stdio.h>
#include <math.h>
#include "PEMstruct.h"	//粒子パラメータを持つ構造体が記述されている

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
	private:int kv;				//動粘性係数
	private:int pe_in_row;
	private:int pe_n;			//粒子数
	private:int wa_n;			//壁の数
	private:double lc;          //代表スケール(セルの大きさ(MD)/壁の位置(DEM))
	private:int step;			//計算ステップ数
	private:double dt;			//計算の時間刻み(s)
	private:PEM *pe;            //粒子パラメータを持つ構造体ポインタ
	private:double camscale;	//カメラ倍率
	private:double *param;		//計算結果の受け渡し用変数ポインタ

	private:Graphics^ gr;		//描画用ポインタ1
	private:Graphics^ gr2;		//描画用ポインタ2
	private:Pen^ blackPen;		//描線用ポインタ（黒）
	private:SolidBrush^ redBrush;	    //塗潰し用ポインタ（赤）
	private:SolidBrush^ blackBrush;	    //塗潰し用ポインタ（黒）
	private:SolidBrush^ blueBrush;	    //塗潰し用ポインタ（青）
	private:SolidBrush^ greenBrush;	    //塗潰し用ポインタ（緑）
	private:array<SolidBrush^>^ atomBrush;	        //各原子の色
	private:BufferedGraphicsContext^ myContext;		//画像バッファ
	private:BufferedGraphics^ myBuffer;				//画像バッファ
	
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
		/// コンストラクタ
		/// クラス生成時に最初に呼ばれる
		/// </summary>
		Form1()
		{
			InitializeComponent();
			//dt=1e-15;                 //計算の時間刻み(s)(MD)
			dt=0.005;					//計算の時間刻み(s)(DEM)
			//lc=1e-9;				    //MDの代表サイズ
			lc=LC;					    //DEMの代表サイズ
			kv=0.1;
			pe_n=5000;					//粒子数の設定
			wa_n=4;						//壁の数の設定
			param = new double[10];		//計算結果などの格納用のメモリ確保
			pe = new PEM[pe_n+wa_n];	//粒子構造体のメモリ確保
			for(int i=0; i<pe_n+wa_n ;i++){
				pe[i].en = new double[pe_n+wa_n];
				pe[i].es = new double[pe_n+wa_n];
				pe[i].coefficient = new double[pe_n+wa_n];
			}
			Config(rows, pe_in_row, wa_n, pe, kv); //粒子の初期値設定(main_MD.cpp参照)

			atomBrush = gcnew array<SolidBrush^>(pe_n);     //原子の色
			redBrush = gcnew SolidBrush( Color::Red );		//塗潰し（赤）
			blackBrush = gcnew SolidBrush( Color::Black );	//塗潰し（黒）
			blueBrush = gcnew SolidBrush( Color::Blue );	//塗潰し（青）
			greenBrush = gcnew SolidBrush( Color::Green );	//塗潰し（青）
			blackPen = gcnew Pen(System::Drawing::Color::Black, 3);	//描線（黒）
			gr = pictureBox1->CreateGraphics();						//描画
		}

	protected:
		/// <summary>
		/// デストラクタ
		/// 使用中のリソースをすべてクリーンアップする
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
		/// デザイナ サポートに必要なメソッドです。このメソッドの内容を
		/// コード エディタで変更しないでください。
		/// 左上の - を押すことでコンパクトにできます。
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

			//プログラムが起動してウィンドウが描画される直前に呼ばれる関数
			//初期設定などしておくとよい
	private: System::Void Form1_Load(System::Object^  sender, System::EventArgs^  e) {
				 //初期値
				 timer1->Interval = 1;     // 何ﾐﾘ秒毎にｲﾍﾞﾝﾄ発生させるか指定
				 timer1->Enabled = false;  // ﾀｲﾏｰを停止状態で初期化
				 step=0;
				 //カメラスケール(=Formの幅÷代表サイズ)
				 camscale=(double)pictureBox1->Width/lc;
				 ScaleLabel->Text = Convert::ToString(camscale);
			 }

			 //STARTボタンがクリックされたとき以下を実行
	private: System::Void StartButton_Click(System::Object^  sender, System::EventArgs^  e) {	
				 label1->Text=trackBar1->Value.ToString("F8");
				 kv = 10^(-1*(trackBar1->Value));
				 Config(rows, pe_in_row, wa_n, pe, kv);
				 timer1->Start();


			 }

			 //STOPボタンをクリックされたとき以下を実行
	private: System::Void StopButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 timer1->Stop();
			 }

			 //スケール拡大
	private: System::Void BiggerButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 camscale *= 1.5;
				 ScaleLabel->Text = Convert::ToString(camscale);
			 }

			 //スケール縮小
	private: System::Void SmallerButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 camscale /= 1.5;
				 ScaleLabel->Text = Convert::ToString(camscale);
			 }

			 //初期条件の再読み込み(RESET)
	private: System::Void ResetButton_Click(System::Object^  sender, System::EventArgs^  e) {
				 for(int i=0; i<pe_n+wa_n ;i++){
					 pe[i].en = new double[pe_n+wa_n];
					 pe[i].es = new double[pe_n+wa_n];
				 }
				 Config(rows, pe_in_row, wa_n, pe, kv); //粒子の初期値設定(PEM.cpp参照)

				 timer1->Interval = 1;     // 何ﾐﾘ秒毎にｲﾍﾞﾝﾄ発生させるか指定
				 timer1->Enabled = false;  // ﾀｲﾏｰを停止状態で初期化
				 step=0;
				 //カメラスケール(=Formの幅÷代表サイズ)
				 camscale=(double)pictureBox1->Width/lc;
				 ScaleLabel->Text = Convert::ToString(camscale);
			 }

			 //タイマーのTickで計算を回す
	private: System::Void timer1_Tick(System::Object^  sender, System::EventArgs^  e) {

				 //粒子法(PEM)のメインプログラム
				 //1Tickで10回計算する
				 for(int i = 0; i < ITERATE_NUM; i++){
					 CalcStep(pe_n, wa_n, pe, dt, lc, param);
					 step++;
				 }

				 //paramに書き込まれたエネルギーなどの情報をForm1に表示


				 /// <summary>
				 /// 以下、画面への描画コード
				 /// </summary>

				 myContext = gcnew BufferedGraphicsContext();	
				 myBuffer = myContext->Allocate(pictureBox1->CreateGraphics(),pictureBox1->DisplayRectangle);
				 myBuffer->Graphics->Clear(Color::White);
				 double camx = (double)pictureBox1->Width;	//pictureBoxの幅
				 double camy = (double)pictureBox1->Height;	//pictureBoxの高さ

				 //粒子の場合
				 for(int ni=0; ni < pe_n+wa_n; ni++){
						 // 表示粒子サイズ(幅)
						 float radius = (float)(camscale*pe[ni].r);
						 // 描画する粒子の位置の指定
						 float xd = (float)(camx/2+camscale*pe[ni].x);
						 float yd = (float)(camy/2-camscale*pe[ni].y);

					 if(pe[ni].type==0){
						 //myBuffer->Graphics->DrawEllipse(blackPen,xd-radius,yd-radius,2*radius,2*radius);
						 myBuffer->Graphics->FillEllipse(blueBrush,xd-radius,yd-radius,2*radius,2*radius);
					 }else{
						 
						 myBuffer->Graphics->FillEllipse(blackBrush,xd-radius,yd-radius,2*radius,2*radius);
					 }
				 }
				 
				 //壁の場合
				 for(int ni=pe_n; ni < pe_n+wa_n; ni++){
					 //壁がないことに設定されていたら飛ばす
					 if(pe[ni].exist == false) continue;

					 //直線の式
					 //la*x + lb*y + lc = 0
					 double la =  sin(pe[ni].phi);
					 double lb = -cos(pe[ni].phi);
					 double lc = -(la*pe[ni].x + lb*pe[ni].y);
					 /*
					 //線の両端のシミュレーション上の座標
					 double x1, y1, x2, y2;

					 if(fabs(lb) < 0.01){ //y軸にだいたい並行な場合
						 y1 = -camx/2.0 / camscale;
						 y2 = camx/2.0 / camscale;

						 x1 = (-lc-lb*y1)/la;
						 x2 = (-lc-lb*y2)/la;
					 }else{
						 x1 = -camx/2.0 / camscale; //x方向両端点の元座標での値
						 x2 = camx/2.0 / camscale;

						 y1 = (-lc-la*x1)/lb;
						 y2 = (-lc-la*x2)/lb;
					 }
					 
					 //カメラ座標にする
					 float camx1 = (float)(camx/2.0 + camscale*x1);
					 float camy1 = (float)(camx/2.0 + camscale*y1);
					 float camx2 = (float)(camx/2.0 + camscale*x2);
					 float camy2 = (float)(camx/2.0 + camscale*y2);

					 //幅
					 float camwidth = (float)(camscale * 2.0 * pe[ni].r);

					 //ペンに幅を設定
					 //C++の文法にない書き方。あまり気にしなくてよい
					 Pen^ linePen = gcnew Pen(System::Drawing::Color::Black, camwidth);
					 //線の描画
					 myBuffer->Graphics->DrawLine(linePen, camx1, camy1, camx2, camy2);
					 */
				 }
				 
				 // 時間の表示
				 double tt = step*dt;
				 TimeLabel->Text=Convert::ToString(tt);
				 myBuffer->Render(gr);
			 }
	private: System::Void trackBar1_Scroll(System::Object^  sender, System::EventArgs^  e) {
				 
			 }
};
}

