#include<stdio.h>
#include<math.h>

double multi(int m,int n,double A[m][n],double w[m],double y[m]){
 int i,j;
 double l[m];
for(i=0;i<m;i++){
       l[i]=0;
      for( j=0;j<m;j++){
      l[i]+=A[i][j]*w[j];}
y[i]=l[i];
}}
double conjugate(int m,int n,double A[m][n],double b[m] ,double x[m],double epsilon,int *count_add,int *count_multi,double resi[100],int*it){
 int i,j,k,l,itr;
double x_old[m],x_new[m],r_old[m],r_new[m],s_old[m],s_new[m],p,q,t[m],u,v,alpha,beta,r_norm,norm,y[m];
       for(i=0;i<m;i++){
           x_old[i]=0;}
multi(m,n,A,x_old,y);
for(i=0;i<m;i++){
    r_old[i]=b[i]-y[i];
     s_old[i]=r_old[i];
     *count_add=*count_add+1;}
itr=1;
*it=1;
do{
multi(m,n,A,s_old,y);
p=0;
q=0;
for(i=0;i<m;i++){
q+=(s_old[i]*y[i]);
p+=(r_old[i]*r_old[i]);
*count_add=*count_add+2;
*count_multi=*count_multi+2;}
alpha=(p/q) ;
*count_multi=*count_multi+1;
for(i=0;i<m;i++){
x_new[i]=x_old[i]+alpha*s_old[i];
   *count_add=*count_add+1;
    *count_multi=*count_multi+1;}
multi(m,n,A,s_old,y);
 for(i=0;i<m;i++){
    r_new[i]=r_old[i]-(alpha*y[i]);
     *count_add=*count_add+1;
    *count_multi=*count_multi+1;}
u=0;v=0;
for(i=0;i<m;i++){
v+=(r_old[i]*r_old[i]);
u+=(r_new[i]*r_new[i]);
*count_add=*count_add+2;
*count_multi=*count_multi+2;}
beta=(u/v) ;
for(i=0;i<m;i++){
s_new[i]=r_new[i]+beta*s_old[i];
*count_add=*count_add+1;
*count_multi=*count_multi+1;}
r_norm=0.0;
for(i=0;i<m;i++){
r_norm+=(r_new[i]*r_new[i]);
*count_add=*count_add+1;
*count_multi=*count_multi+1;}
for(i=0;i<m;i++){
x_old[i]=x_new[i];
r_old[i]=r_new[i];
s_old[i]=s_new[i];}
norm=sqrt(r_norm);
resi[itr]=norm;
*it=*it+1;
itr++;
}while(norm>epsilon);
for(l=0;l<m;l++){
        x[l]=x_new[l];}}

int main(){
 int i,j,count_add,count_multi,it;
FILE *f[2],*f_out,*f_q2,*f_outplot,*f_out_1b,*f_q2_2b;
int m[2],n[2],p,q,r,s,l;
char c;
for(i=0; i<2; i++){
        if(i==0){
f[i]=fopen("input_1a.txt","r");
}
if(i==1){
f[i]=fopen("input_1b.txt","r");
}
         m[i]=0;n[i]=1;

while((c=fgetc(f[i]))!=EOF){
  if(c=='\t'&&m[i]==0){
    n[i]++;}
  else if(c=='\n'){
  m[i]++;}
}}
p=m[0]; q=n[0];
r=m[1]; s=n[1];
fclose(f[0]);
fclose(f[1]);
double A[p][q],C[r][s],b[p],d[r],x1[p],x2[r],resi[100];
f[0]=fopen("input_1a.txt","r");
for( i=0; i<p; i++){
    for( j=0; j<q ; j++){
        fscanf(f[0],"%lf",&A[i][j]);}}
fclose(f[0]);
f[1]=fopen("input_1b.txt","r");
for( i=0; i<r; i++){
    for( j=0; j<s ; j++){
        fscanf(f[1],"%lf",&C[i][j]);}}
fclose(f[1]);
printf("\nAugumented matrix of Q_1_a\n");

f_out=fopen("Quetion_1_output.txt","w");
f_outplot=fopen("residual_vs_itr.txt","w");
f_q2=fopen("Question_2_output.txt","w");
f_out_1b=fopen("Quetion_b_output.txt","w");
f_q2_2b=fopen("Question_bb_output.txt","w");
fprintf(f_out,"\nAugumented matrix of Q_1_a\n");
for( i=0;i<p;i++){
       for( j=0;j<q;j++){
       printf("%lf\t",A[i][j]);
       fprintf(f_out,"%f\t",A[i][j]);
        fprintf(f_q2,"%f\t",A[i][j]);
        }
       printf("\n");
       fprintf(f_out,"\n");
    fprintf(f_q2,"\n");}
fprintf(f_out,"\n");
fprintf(f_q2,"\n");
for( i=0;i<p;i++){
    j=q-1;
    b[i]=A[i][j];
}
for( i=0;i<r;i++){
    j=s-1;
    d[i]=C[i][j];
}
double epsilon;
fprintf(f_outplot,"\n\n**r_values**\n***residual values***\n\n");
for(epsilon=pow(10,-3);epsilon>=pow(10,-6);epsilon=epsilon/10.0){
count_add=0;  count_multi=0;
conjugate(p,q,A,b,x1,epsilon,&count_add,&count_multi,resi,&it);
 fprintf(f_outplot,"epsilon= %lf\n",epsilon);
for(i=1;i<=(it-1);i++){
    //printf("\nr value :%lf",resi[i]);
  fprintf(f_outplot,"\n%d\t%e",i,resi[i]);
}
fprintf(f_outplot,"\n\n");
printf("\nepsilon = %lf\t:\n",epsilon);
fprintf(f_out,"epsilon = %lf\t:\n",epsilon);
fprintf(f_out,"\nsolution is\n");
printf("\nsolution is\n");
fprintf(f_q2,"epsilon = %lf\t:\n",epsilon);
fprintf(f_q2,"\nsolution is\n");
for(l=0;l<p;l++){
printf("%0.9lf\n",x1[l]);
fprintf(f_out,"%0.6lf\n",x1[l]);
fprintf(f_q2,"%0.6lf\n",x1[l]);
}
fprintf(f_out,"\n\n\n");
printf("number of addition and subtractions are  %d\n",count_add);
printf("number of multiplications and divisions are  %d\n\n",count_multi);
fprintf(f_q2,"number of addition and subtractions are  %d\n",count_add);
fprintf(f_q2,"number of multiplications and divisions are  %d\n\n",count_multi);
}
printf("\n\n");
fclose(f_out);
fclose(f_outplot);
fclose(f_q2);
printf("\nAugumented matrix of Q_1_b\n");
fprintf(f_out,"\nAugumented matrix of Q_1_b\n");
for( i=0;i<r;i++){
       for( j=0;j<s;j++){
       printf("%lf\t",C[i][j]);
       fprintf(f_out,"%f\t",C[i][j]);
        fprintf(f_q2,"%f\t",C[i][j]);
        }
       printf("\n");
       fprintf(f_out,"\n");
     fprintf(f_q2,"\n");
}
fprintf(f_out,f_q2,"\n");
for(epsilon=pow(10,-3);epsilon>=pow(10,-6);epsilon=epsilon/10.0){
         count_add=0;
         count_multi=0;
conjugate(r,s,C,d,x2,epsilon,&count_add,&count_multi,resi,&it);
printf("epsilon = %lf\t:\n",epsilon);
fprintf(f_out_1b,"epsilon = %lf\t:\n",epsilon);
fprintf(f_q2_2b,"epsilon = %lf\t:\n",epsilon);
for(l=0;l<r;l++){
printf("%0.6lf\n",x2[l]);
fprintf(f_out_1b,"%0.6lf\n",x2[l]);
fprintf(f_q2_2b,"%0.6lf\n",x2[l]);
}
printf("\n");
 fprintf(f_out_1b,"\n");
 printf("number of addition and subtractions are  %d\n",count_add);
  printf("number of multiplications and divisions are  %d\n\n",count_multi);
  fprintf(f_q2_2b,"number of addition and subtractions are  %d\n",count_add);
fprintf(f_q2_2b,"number of multiplications and divisions are  %d\n\n",count_multi);
}
}
