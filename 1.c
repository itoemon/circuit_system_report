#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define N 1000//サンプリング数
#define Wc 5//カットオフ周波数

/* LPF関数(1次，バタワース) */
void LPF(double *ReH, double *ImH){
  int k;
  double W;
  for(k=0;k<N;k++){
    W = (double)k/Wc;
    ReH[k] = ReH[k] / (1.0+W*W);
    ImH[k] = ImH[k] * (W/(1.0+W*W));
  }
}

int main(int argc, char *argv[]){
  double input[N];//入力信号
  double output1[N],output2[N];//出力信号
  double Re[N],Im[N];//DFT後の実部成分と虚部成分
  double A;//LPF伝達関数計算用
  double ReH[N],ImH[N];//LPF伝達関数の実部成分と虚部成分
  double Z1,Z2;//遅延器出力信号
  double Qout;//1bit量子化出力
  double tmp,tmp2;//格納用
  double Spec;//振幅スペクトル出力用
  FILE *fp;
  int i,j;
  double W = 2.0*M_PI/N;

  //LPF伝達関数の振幅スペクトル計算
  for(i=0;i<N;i++){ 
    A = (double)i/Wc;
    ReH[i] = 1.0/(1.0+A*A);
    ImH[i] = (A/(1.0+A*A));
  }
  //LPF伝達関数の振幅スペクトル出力
  fp = fopen("H_spec.txt", "wb");
  for(i=0; i<N; i++){
    Spec = sqrt(ReH[i]*ReH[i]+ImH[i]*ImH[i]);
    fprintf(fp, "%d %f\n", i, Spec);
  }

/*** 入力信号のサンプリング ***/
  //初期化
  for(i=0; i<N; i++){
    input[i] = 0.0;
    Re[i] = 0.0;
    Im[i] = 0.0;
    ReH[i] = 0.0;
    ImH[i] = 0.0;
  }
  //オーバーサンプリング
  for(i=0; i<N; i++){
    input[i] = sin(W*i);
  }
  //入力信号の時間波形出力
  fp = fopen("Input_time.txt", "wb");
  for(i=0; i<N; i++){
    fprintf(fp, "%d %f\n", i, input[i]);
  }
  //入力信号を実部と虚部に分けてDFT
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      Re[i] += input[j]*cos(W*i*j);
      Im[i] += -input[j]*sin(W*i*j);
    }
  }
  //入力信号の振幅スペクトル出力
  fp = fopen("Input_spec.txt", "wb");
  for(i=0; i<N; i++){
    Spec = sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
    fprintf(fp, "%d %f\n", i, Spec);
  }

/*** １次のΔΣ変調 ***/
  //初期化
  for(i=0; i<N; i++){
    output1[i] = 0.0;
    Re[i] = 0.0;
    Im[i] = 0.0;
  }
  tmp = 0.0;
  //１次のΔΣ変調計算部分
  for(i=0; i<N; i++){
    //１つ前の出力を足す
    if(i!=0){
      tmp = input[i] + Z1;
    }
    //1bit量子化
    if(tmp>=0.0){
      Qout = 1.0;
    }else if(tmp<0.0){
      Qout = -1.0;
    }
    //出力とフィードバック処理
    output1[i] = Qout;
    Z1 = tmp - output1[i];
  }
  //変調信号の時間波形出力
  fp = fopen("1stout_time.txt", "wb");
  for(i=0; i<N; i++){
    fprintf(fp, "%d %f\n", i, output1[i]);
  }
  //変調信号を実部と虚部に分けてDFT
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      Re[i] += output1[j]*cos(W*i*j);
      Im[i] += -output1[j]*sin(W*i*j);
    }
  }
  //変調信号の振幅スペクトル出力
  fp = fopen("1stout_spec.txt", "wb");
  for(i=0; i<N; i++){
    Spec = sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
    fprintf(fp, "%d %f\n", i, Spec);
  }
  //LPFを適用
  LPF(Re, Im);
  //LPFを適用した後の振幅スペクトル出力
  fp = fopen("1stout_LPF_spec.txt", "wb");
  for(i=0; i<N; i++){
    Spec = sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
    fprintf(fp, "%d %f\n", i, Spec);
  }
  //IDFT
  for(i=0;i<N;i++){
    output1[i] = 0.0;//初期化
    for(j=0;j<N;j++){
      output1[i] += Re[j]*cos(2*M_PI*j*i/N) - Im[j]*sin(2*M_PI*j*i/N);
    }
   output1[i] /= N;
  }
  //LPFした後の時間波形出力
  fp = fopen("1stout_LPF_time.txt", "wb");
  for(i=0; i<N; i++){
    fprintf(fp, "%d %f\n", i, output1[i]);
  }
  
/*** 2次のΔΣ変調 ***/
  //初期化
  for(i=0; i<N; i++){
    output2[i] = 0.0;
    Re[i] = 0.0;
    Im[i] = 0.0;
  }
  tmp = 0.0;
  tmp2 = 0.0;
  //2次のΔΣ変調計算部分
  for(i=0; i<N; i++){
    //2つ前の出力を足す
    if(i>=2){
      tmp = input[i] - Z2;
    }
    //1つ前の出力を足す
    if(i>=1){
      tmp2 = tmp + 2.0*Z1;
    }
    //1bit量子化
    if(tmp2>=0.0){
      Qout = 1.0;
    }else if(tmp2<0.0){
      Qout = -1.0;
    }
    //出力とフィードバック処理
    output2[i] = Qout;
    Z1 = tmp2 - output2[i];
    Z2 = Z1;
  }
  //変調信号の時間波形出力
  fp = fopen("2ndout_time.txt", "wb");
  for(i=0; i<N; i++){
    fprintf(fp, "%d %f\n", i, output2[i]);
  }
  //変調信号を実部と虚部に分けてDFT
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      Re[i] += output2[j]*cos(W*i*j);
      Im[i] += -output2[j]*sin(W*i*j);
    }
  }
  //変調信号の振幅スペクトル出力
  fp = fopen("2ndout_spec.txt", "wb");
  for(i=0; i<N; i++){
    Spec = sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
    fprintf(fp, "%d %f\n", i, Spec);
  }
  //LPFを適用
  LPF(Re, Im);
  //LPFした後の振幅スペクトル出力
  fp = fopen("2ndout_LPF_spec.txt", "wb");
  for(i=0; i<N; i++){
    Spec = sqrt(Re[i]*Re[i]+Im[i]*Im[i]);
    fprintf(fp, "%d %f\n", i, Spec);
  }
  //IDFT
  for(i=0;i<N;i++){
    output2[i] = 0.0;//初期化
    for(j=0;j<N;j++){
      output2[i] += Re[j]*cos(2*M_PI*j*i/N) - Im[j]*sin(2*M_PI*j*i/N);
    }
   output2[i] /= N;
  }
  //LPFした後の時間波形出力
  fp = fopen("2ndout_LPF_time.txt", "wb");
  for(i=0; i<N; i++){
    fprintf(fp, "%d %f\n", i, output2[i]);
  }
fclose(fp);
}
