#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MEM 100
#define NMAX 100
#define PI 3.141592

double EPS_MIN = 0.000000001;
double SEA = 100;
double SEA_H = 0.1;
double OVER_Y = 10000;
double OVER_X = 0.001;
double H = 0.001;
double ANSMAX = 100;
FILE *gs;

double x,xo,error,fx,fxd,xph,xmh,fxp,fxm,ini_count;
double fx_ch[MEM],fxd_ch[MEM],x_ch[MEM],fxdo;
int i,j,p,e,count,match,ansn,y_over,x_over,mode[MEM];

void rivers(){
    if(ini_count<=0){
        ini_count-=SEA_H;
    }else{
        ini_count+=SEA_H;
    }
    ini_count=-ini_count;
    xo=ini_count;
}
void ans_area(){
    count=0;
}
void inii(){
    ini_count=0;
    x=ini_count;
    mode[1]=1;
    mode[4]=0;
    count=0;
}
double der(double fxp,double fxm){
    return((fxp-fxm)/(2*H));
}
double newton_raphson(double fx,double fxd){
    return(x-(fx/fxd));
}
void diver(){
    if(count==1){
        if(mode[3]==1){
            fx_ch[0]=fx;
            fx_ch[1]=fabs(fx_ch[0]*x);
            fx_ch[3]=fabs(fx_ch[0]-fx_ch[2]);
            if(fx_ch[1] > OVER_Y){
                y_over++;
            }else{
                y_over=0;
            }
            if(fx_ch[3] < OVER_X){
                x_over++;
            }else if(isnan(fx_ch[3])==0){
                x_over=0;
            }
            if(y_over==10 || x_over==10){
                mode[4]=1;
            }
        }else if(mode[3]==2){
            if(ini_count>=SEA){
                mode[4]=1;
            }
        }else{}
        if(ansn>=ANSMAX){
            mode[4]=1;
        }
        if(mode[4]==1){
            printf("\n>> 計算が終了しました。\n\n");
            printf(">> 解の個数 : %d\n\n",ansn);
            for(i=0;i<ansn;i++){
                printf(">> x[%d] = %3.8f\n",i+1,x_ch[i]);
            }
            mode[1]=0;
            for(i=0;i<MEM;i++){
                x_ch[i]=0;
            }
            y_over=0;
            x_over=0;
            ansn=0;
            mode[1]=0;
            mode[4]=0;
            printf("\n何かキーを入力してください。");
            getchar();
            getchar();
        }
        if(isnan(fx_ch[0])){
        }else{
            fx_ch[2]=fx_ch[0];
        }
    }
}
void eps(){
    error = fabs(xo-x);
    if(error<=EPS_MIN){
        pclose(gs);
        if(mode[3]==1){
            x_ch[ansn]=x;
            for(i=0;i<=ansn;i++){
                if(fabs(x_ch[ansn]-x_ch[i])<EPS_MIN){
                    match++;
                }
            }
            ansn++;
            count=0;
            rivers();
            if(match>1){
                ansn-=1;
            }
        }else if(mode[3]==2){
            x_ch[ansn]=x;
            for(i=0;i<=ansn;i++){
                if(fabs(x_ch[ansn]-x_ch[i])<EPS_MIN){
                    match++;
                }
            }
            ansn++;
            count=0;
            rivers();
            if(match>1){
                ansn-=1;
            }
        }else if(mode[3]==3){
            x_ch[0]=x;
            ansn=1;
            mode[4]=1;
            printf("\n>> %d回で収束しました。",count-1);
            count=0;
        }
    }else{
        mode[4]=0;
        if(mode[3]==1){
            if(count>NMAX){
                count=0;
                rivers();
            }
        }else if(mode[3]==2){
            if(count>NMAX){
                count=0;
                rivers();
            }
        }else if(mode[3]==3){
            if(count>NMAX){
                count=0;
                printf("(警告)この初期値では収束しません。\n");
                mode[4]=0;
                mode[1]=0;
                printf("\n何かキーを入力してください。");
                getchar();
                getchar();
            }
        }
    }
    fxdo=fxd;
    match=0;
    x = xo;
}
double func(double x){
    double y;
    if(mode[0]==1){
        y=x+cos(x);
    }else if(mode[0]==2){
        y=pow(x,2)-4;
    }else if(mode[0]==3){
        y=4*pow(x,3)-2*pow(x,2)-6*x+3;
    }else if(mode[0]==4){
        y=pow(x,3)-2*pow(x,2)-x+2;
    }else if(mode[0]==5){
        y=pow(x,4)-7*pow(x,2)+12;
    }else if(mode[0]==6){
        y=(x*x*x*x)-(13*x*x)+36;
        y=pow(x,4)-13*pow(x,2)+36;
    }else if(mode[0]==7){
        y=exp(x)+x;
    }else if(mode[0]==8){
        y=pow(x,4)-pow(x,3)-pow(x,2)-2*x+1;
    }else if(mode[0]==9){
        y=(x-cos(x))/(x+5+log(x));
    }else if(mode[0]==10){
        y=pow(x,2)-2*x;
    }else if(mode[0]==11){
        y=pow(x,5)-5*pow(x,4)+8*pow(x,3)+32*pow(x,2)-149*x+113;
    }else if(mode[0]==12){
        y=pow(x,6)-21*pow(x,5)+175*pow(x,4)-735*pow(x,3)+1624*pow(x,2)-1764*x+720;
    }else{}
    return(y);
}
int main(void){
    mode[3]=3;
    while(1){
        inii();
        puts("========================= 情報表示 ============================");
        puts("【設定情報】");
        if(mode[3]==1){
            puts("   現在の設定 : 区間解析(すべての解を探索)");
            puts("   プロットモード : ");
            if(mode[6]==1){
                puts("ON");
            }else if(mode[6]==0){
                puts("OFF");
            }
            puts("【各パラメータ情報】");
            printf("   収束解の最小許容誤差 : %3.13lf\n",EPS_MIN);
            printf("   Y軸方向への発散判定とする値 : %3.13lf\n",OVER_Y);
            printf("   X軸方向への発散判定とする値 : %3.13lf\n",OVER_X);
            printf("   区間探索精度(小さいほど精度は高い) : %3.13lf\n",SEA_H);
            printf("   数値微分精度(小さいほど精度は高い) : %3.13lf\n",H);
        }else if(mode[3]==2){
            puts("   現在の設定 : 手動区間指定");
            puts("   プロットモード : ");
            if(mode[6]==0){
                puts("ON");
            }else if(mode[6]==1){
                puts("OFF");
            }
            puts("【各パラメータ情報】");
            printf("   収束解の最小許容誤差 : %3.13lf\n",EPS_MIN);
            printf("   区間指定幅 : %3.13lf\n",SEA);
            printf("   区間探索精度(小さいほど精度は高い) : %3.13lf\n",SEA_H);
            printf("   数値微分精度(小さいほど精度は高い) : %3.13lf\n",H);
        }else if(mode[3]==3){
            puts("   現在の設定 : 初期値指定");
            printf("   プロットモード : ");
            if(mode[6]==0){
                puts("ON");
            }else if(mode[6]==1){
                puts("OFF");
            }
            puts("【各パラメータ情報】");
            printf("   収束解の最小許容誤差 : %3.13lf\n",EPS_MIN);
        }
        puts("========================== 式一覧 ===========================");
        puts("【1】  (x)+cos(x)");
        puts("【2】  (x^2)-4");
        puts("【3】  (4x^3)-(2x^2)-(6x)+3");
        puts("【4】  (x^3)-(2x^2)-(x)+2");
        puts("【5】  (x^4)-(7x^2)+12");
        puts("【6】  (x^4)-(13x^2)+36");
        puts("【7】  (e^x)+x");
        puts("【8】  (x^4)-(x^3)-(x^2)-(2x)+1");
        puts("【9】  (x-cos(x))/(x+5+log(x))");
        puts("【10】 (x^2)-(2x)");
        puts("【11】 (x^5)-(5x^4)+(8x^3)+(32x^2)-(149x)+113");
        puts("【12】 (x^6)-(21x^5)+(175x^4)-(735x^3)+(1624x^2)-1764x+720");
        puts("【20】 (解析法選択)");
        puts("【21】 (パラメータ設定)");
        puts("【22】 (終了)");
        puts("【23】 (プロットモード)");
        puts("=============================================================");
        printf("番号を入力してください => ");
        scanf("%d",&mode[0]);
        puts("=============================================================");
        if(mode[0]==22){
            puts("See you!!");
            break;
        }else if(mode[0]==23){
            if(mode[6]==0){
                mode[6]=1;
            }else if(mode[6]==1){
                mode[6]=0;
            }
            mode[1]=0;
        }else if(mode[0]==21){
            puts("パラメータ設定");
            puts("パラメータ設定を行います。");
            printf("【1】収束解の最小許容誤差 : %3.13lf\n",EPS_MIN);
            printf("【2】Y軸方向への発散判定とする値 : %3.13lf\n",OVER_Y);
            printf("【3】X軸方向への発散判定とする値 : %3.13lf\n",OVER_X);
            printf("【4】区間探索精度(小さいほど精度は高い) : %3.13lf\n",SEA_H);
            printf("【5】数値微分精度(小さいほど精度は高い) : %3.13lf\n",H);
            printf("【6】区間指定幅 : %3.13lf\n",SEA);
            puts("=============================================================");
            printf("番号を入力してください => ");
            scanf("%d",&mode[5]);
            if(mode[5]==1){
                printf("収束解の最小許容誤差 = ");
                scanf("%lf",&EPS_MIN);
            }else if(mode[5]==2){
                printf("Y軸方向への発散判定とする値 = ");
                scanf("%lf",&OVER_Y);
            }else if(mode[5]==3){
                printf("X軸方向への発散判定とする値 = ");
                scanf("%lf",&OVER_X);
            }else if(mode[5]==4){
                printf("区間探索精度(小さいほど精度は高い) = ");
                scanf("%lf",&SEA_H);
            }else if(mode[5]==5){
                printf("数値微分精度(小さいほど精度は高い) = ");
                scanf("%lf",&H);
            }else if(mode[5]==6){
                printf("区間指定幅 = ");
                scanf("%lf",&SEA);
            }else{}
            mode[1]=0;
        }else if(mode[0]==20){
            puts("【1】自動解析モード(すべての解を探索)");
            puts("【2】区間指定モード");
            puts("【3】初期値指定モード");
            puts("=============================================================");
            printf("番号を入力してください => ");
            scanf("%d",&mode[3]);
            mode[1]=0;
        }else{
            if(mode[3]==3){
                printf("初期値を入力してください => ");
                scanf(" %lf",&x);
            }
        }
        gs = popen("gnuplot -persist","w");
        fprintf(gs, "plot '-' with linespoints pointtype 2\n");
        while(mode[1]==1){
            fx=func(x);
            diver();
            fxp=func(x+H);
            fxm=func(x-H);
            fxd=der(fxp,fxm);
            //printf("fxd = %f...fx = %f\n",fxd,fx);
            xo=newton_raphson(fx,fxd);
            eps();
            if(mode[6]==0 && mode[3]==1 && mode[4]==0){
                fprintf(gs, "%lf\n", error);
            }
            if(mode[6]==0 && mode[3]==3 && mode[4]==0){
                fprintf(gs, "%lf\n", error);
                printf("%d回目 = %lf\n",count,error);
            }
            count++;
        }
        pclose(gs);
    }
}
