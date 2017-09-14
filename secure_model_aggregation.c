#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <obliv.h>
#include "secure_model_aggregation.h"
#include <time.h>

#include"../common/util.h"

#define NUM 1

int main(int argc,char *argv[])
{
  char c[30], file1[50],file2[50], file3[50], file4[50];
  char str[15];

  int i, j;
  FILE *fptr1, *fptr2, *fptr3, *fptr4;
  ProtocolDesc pd;
  protocolIO io;
  clock_t start, end;
  double cpu_time_used;

  strcpy(file1,"Inputs/");
  strcat(file1, str);
  file1[strlen(file1)] = '\0';
  strcat(file1, "_beta1.txt");
  
  strcpy(file2,"Inputs/");
  strcat(file2, str);
  file2[strlen(file2)] = '\0';
  strcat(file2, "_beta2.txt");

  strcpy(file3,"Results/20 parties/");
  strcat(file3, str);
  file3[strlen(file3)] = '\0';
  strcat(file3, ".txt");

  strcpy(file4,"Results/20 parties/time");
  strcat(file4, str);
  file4[strlen(file4)] = '\0';
  strcat(file4, ".txt");

  if(argv[1][0]=='1')
  {
    fptr1=fopen(file1,"r");
    for(i = 0; i < M; i++)
    {
      for(j = 0; j < D; j++)
      {
        fscanf(fptr1,"%s", c);
        io.beta1[i][j] = atoi(c);
      }
    }
    fclose(fptr1);
  }
  if(argv[1][0]=='2')
  {
    fptr2=fopen(file2,"r");
    for(i = 0; i < M; i++)
    {
      for(j = 0; j < D; j++)
      {
        fscanf(fptr2,"%s", c);
        io.beta2[i][j] = atoi(c);
      }
    }
    fclose(fptr2);
  }

  start = clock();
  protocolUseStdio(&pd);
  setCurrentParty(&pd, argv[1][0] == '1' ? 1 : 2);
  execYaoProtocol(&pd, aggregate, &io);
  cleanupProtocol(&pd);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  fprintf(stderr,"\nElapsed Time: %f\n",cpu_time_used);

  if(argv[1][0]=='1')
  {
    fptr3=fopen(file3,"w");
    for(i = 0; i < D; i++) 
       fprintf(fptr3,"%.8f ",io.beta_avg[i]*1.0/SCALE);
    fclose(fptr3);
  }
    
  fptr4=fopen(file4,"a");
  fprintf(fptr4,"Elapsed Time: %f\n",cpu_time_used);
  fclose(fptr4);

  return 0;
}
