#include<stdio.h>



double Gauss_Seidel(int m,int n,double a[m][n],double b[m] ,double x_gauss[m],double epsilon ,int *count_add,int *count_multi){

    int i,j,k,l,itr;
    double x[m],er[m],sum,temp,max_er;

     for(i=0;i<m;i++){

     x[i]=0;
}
itr=1;
do{
for(i=0;i<m;i++){
     sum=0 ;
    for(j=0;j<m;j++){
             if(i !=j){
           sum+=a[i][j]*x[j];
           *count_add=*count_add+1;
              *count_multi=*count_multi+1;
        }
}

            temp=x[i];
            x[i]=(1/a[i][i])*(b[i]-sum);
            er[i]=((x[i]-temp)/x[i]) ;
            *count_add=*count_add+2;
            *count_multi=*count_multi+3;


    }
max_er=0;
for(i=0;i<m;i++){
if(fabs(er[i])>max_er){
    max_er=er[i];
    }
}

itr+=1;

}while(max_er>epsilon);

for(i=0;i<m;i++){
    x_gauss[i]=x[i];
}
}



int main(){

   int i,j,count_add,count_multi;
FILE *f[2],*f_out,*f_q3;
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
//printf("%d\t%d\n",p,q);
r=m[1]; s=n[1];
//printf("%d\t%d\n",r,s);
 fclose(f[0]);
 fclose(f[1]);
double A[p][q],C[r][s],b[p],d[r],x1[p],x2[r];
f[0]=fopen("input_1a.txt","r");
for( i=0; i<p; i++){
    for( j=0; j<q ; j++){
        fscanf(f[0],"%lf",&A[i][j]);
    }}
fclose(f[0]);
f[1]=fopen("input_1b.txt","r");
for( i=0; i<r; i++){
    for( j=0; j<s ; j++){
        fscanf(f[1],"%lf",&C[i][j]);
    }}
fclose(f[1]);
for( i=0;i<p;i++){
    j=q-1;
    b[i]=A[i][j];
}
for( i=0;i<r;i++){
    j=s-1;
    d[i]=C[i][j];
}
printf("\nAugumented matrix of Q_2_a\n");
printf(f_out,"\nAugumented matrix of Q_2_a\n");
for( i=0;i<p;i++){
       for( j=0;j<q;j++){
       printf("%lf\t",A[i][j]);
        }
       printf("\n");}
 printf("\n");
 f_out=fopen("Quetion_2_output.txt","w");
 f_q3=fopen("Question_3_output_Gauss.txt","w");
fprintf(f_out,"Gauss Siedel\n\n");

for( i=0; i<p; i++){
    for( j=0; j<q ; j++){
        fprintf(f_out,"%f\t",A[i][j]);
        fprintf(f_q3,"%f\t",A[i][j]);
    }
    fprintf(f_out,"\n");
    fprintf(f_q3,"\n");
    }
fprintf(f_out,"\n");
 fprintf(f_q3,"\n");

double epsilon;
for(epsilon=pow(10,-3);epsilon>=pow(10,-6);epsilon=epsilon/10.0){
        count_add=0;
        count_multi=0;
       Gauss_Seidel(p,q,A,b,x1,epsilon,&count_add,&count_multi);
       printf("\nresidual error = %lf\t:\n",epsilon);
       fprintf(f_out,"\nresidual error = %lf\t:\n",epsilon);
       fprintf(f_q3,"\n\nresidual error = %lf\t:\n",epsilon);
  for(l=0;l<p;l++){
      printf("%0.6lf\n",x1[l]);
    }
   printf("\n");
    for(l=0;l<p;l++){
fprintf(f_out,"%0.6lf\n",x1[l]);
fprintf(f_q3,"%0.6lf\n",x1[l]);
}
fprintf(f_out,"\n\n\n");
printf("number of addition and subtractions are  %d\n",count_add);
printf("number of multiplications and divisions are  %d\n\n\n",count_multi);
fprintf(f_q3,"number of addition and subtractions are  %d\n",count_add);
fprintf(f_q3,"number of multiplications and divisions are  %d\n\n",count_multi);
}
printf("\n");
printf("\nAugumented matrix of Q_2_b\n");
fprintf(f_out,"\nAugumented matrix of Q_2_b\n");
for( i=0;i<r;i++){
       for( j=0;j<s;j++){
       printf("%lf\t",C[i][j]);
        }
       printf("\n");

}
 printf("\n\n");
 for( i=0; i<r; i++){
    for( j=0; j<s; j++){
        fprintf(f_out,"%f\t",C[i][j]);
        fprintf(f_q3,"%f\t",C[i][j]);
    }
    fprintf(f_out,"\n");
    fprintf(f_q3,"\n");
    }
fprintf(f_out,"\n");
fprintf(f_q3,"\n");

for(epsilon=pow(10,-3);epsilon>=pow(10,-6);epsilon=epsilon/10.0){
        count_add=0;
        count_multi=0;
Gauss_Seidel(r,s,C,d,x2,epsilon,&count_add,&count_multi);
printf("\nresidual error = %lf\t:\n",epsilon);
fprintf(f_out,"\nresidual error = %lf\t:\n",epsilon);
fprintf(f_q3,"\n\nresidual error = %lf\t:\n",epsilon);
for(l=0;l<r;l++){
printf("%0.6lf\n",x2[l]);

}
for(l=0;l<r;l++){
fprintf(f_out,"%0.6lf\n",x2[l]);
fprintf(f_q3,"%0.6lf\n",x2[l]);


}
 //printf("\n");
 //fprintf(f_out,"\n");
 printf("number of addition and subtractions are  %d\n",count_add);
printf("number of multiplications and divisions are  %d\n\n",count_multi);
fprintf(f_q3,"number of addition and subtractions are  %d\n",count_add);
fprintf(f_q3,"number of multiplications and divisions are  %d\n\n",count_multi);
}
printf("\n");

}


