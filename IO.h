#include <stdio.h>

#define INPUTDATA(Type, flag)                                                                             \
void InputData_##Type(Type * input, const int size, const char * name_file){\
  FILE * fp = fopen(name_file, "r");                                                                \
  if (fp==NULL){                                                                                    \
    printf("open file: %s failed!", name_file);                                                     \
    exit(-1);                                                                                       \
  }                                                                                                 \
  int cnt=0;                                                                                        \
  while( cnt< size ){                                                                               \
    if( fscanf(fp,#flag, &input[cnt]) >0 ){                                                          \
      cnt++;                                                                                        \
    }                                                                                               \
    else{                                                                                           \
      break;                                                                                        \
    }                                                                                               \
  }                                                                                                 \
  if(cnt<size){                                                                                     \
    printf("from file %s : got %d raw data, which is less than size- %d !\n", name_file, cnt, size); \
  }                                                                                                 \
  fclose(fp);                                                                                       \
}

#define OUTPUTDATA(Type, flag)                                                                              \
void OutputData_##Type(Type * output, const int size, const char * name_file){                              \
  FILE * fp = fopen(name_file, "w");                                                                  \
  if (fp==NULL){                                                                                      \
    printf("open file: %s failed!", name_file);                                                       \
    exit(-1);                                                                                         \
  }                                                                                                   \
  int cnt=0;                                                                                          \
  while( cnt< size ){                                                                                 \
    fprintf(fp,#flag,output[cnt]);                                                                    \
    fprintf(fp," ");                                                                    \
    cnt++;                                                                                            \
  }                                                                                                   \
  if(cnt<size){                                                                                       \
    printf("write %d data into file %s , which is less than size- %d !\n", cnt, name_file, size);      \
  }                                                                                                   \
  fclose(fp);                                                                                         \
}


INPUTDATA(float, %f) 
INPUTDATA(double, %lf) 
INPUTDATA(int, %d) 
OUTPUTDATA(float, %f)
OUTPUTDATA(double, %f)
OUTPUTDATA(int, %d)









