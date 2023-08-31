#include<stdio.h>

double inverse(int n,double B[n][2*n],double A_inverse[n][n]){
double p;int i,j;
for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){
        if(i>j){
                if(B[j][j]==0){
            for(int t=j+1 ; t<n ; t++){
                if(B[t][j]!=0){
                        for(int l=0;l<2*n;l++){
                            float tempvar;
                            tempvar=B[j][l];
                            B[j][l]=B[t][l];
                            B[t][l]=tempvar;}}}
            }
            float val=(B[i][j]/B[j][j]);
            for(int k=0;k<2*n;k++){
                B[i][k]=B[i][k]-val*B[j][k];}}}}
for(int i=0;i<n;i++){
        p=B[i][i];
   for(int j=0;j<2*n;j++){
    B[i][j]=B[i][j]/p;
}}
for(int j=(n-1);j>=0;j--){
         for(int i=n-1;i>=0;i--){
    if(i<j){
    float val=(B[i][j]/B[j][j]);
            for(int k=0;k<2*n;k++){
                B[i][k]=B[i][k]-val*B[j][k];
                  }}}}
       for(int i=0;i<n;i++){
       for(int j=n;j<2*n;j++){
         A_inverse[i][j-n]=B[i][j];
        }}}

double eigen(int n,double A[n][n],double *eigenvalue){
int i,iter=0,j;double x_old[n],y[n],max_y,x[n],er,maxold_y;
for(i=0;i<n;i++){
x_old[i]=1.0;}
maxold_y=0;
er=1;
do{
for(i=0;i<n;i++){
y[i]=0;
for(j=0;j<n;j++){
  y[i]+=A[i][j]*x_old[j];}}
max_y=0;
for(i=0;i<n;i++){
if(fabs(y[i])>max_y){
    max_y=fabs(y[i]);}}
for(i=0;i<n;i++){
x[i]=(y[i]/max_y);
x_old[i]=x[i];}
er=fabs((max_y-maxold_y));
maxold_y=max_y;
iter++;
}while(er>1e-9);
*eigenvalue=max_y;
}
int main(){
int n;
n=4;FILE *f_out;
double A_inverse[n][n],eigenvalue;
double A[4][4]={{5, 4, 1, 1},{4, 5, 1, 1},{1, 1, 4, 2},{1, 1, 2, 4}};
double B[4][8]={{5, 4, 1, 1,1,0,0,0},{4, 5, 1, 1,0,1,0,0},{1, 1, 4, 2,0,0,1,0},{1, 1, 2, 4,0,0,0,1}};
f_out=fopen("output.txt","w");
printf("\nMatrix A is\n");
fprintf(f_out,"\nMatrix A is\n");
for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        printf("%f\t",A[i][j]);
        fprintf(f_out,"%f\t",A[i][j]);}
    printf("\n");
    fprintf(f_out,"\n");
    }
printf("\n");
eigen(n,A,&eigenvalue);
printf("\n maximum eigenvalue is : %0.6lf\n\n",eigenvalue);
fprintf(f_out,"\n maximum eigenvalue is : %0.6lf\n",eigenvalue);
inverse(n,B,A_inverse);
eigen(n,A_inverse,&eigenvalue);
printf(" minimum eigenvalue is : %0.6lf\n\n",(1/eigenvalue));
fprintf(f_out," minimum eigenvalue is : %0.6lf\n\n",(1/eigenvalue));
}

