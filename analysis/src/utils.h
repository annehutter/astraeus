#ifndef UTILS_H
#define UTILS_H

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef MAXLENGTH
#define MAXLENGTH 1024
#endif

/* ------------------------------------------------------------*/
/* RANDOM NUMBERS                                              */
/* ------------------------------------------------------------*/
double randdouble();

/* ------------------------------------------------------------*/
/* STRING ADDITION                                             */
/* ------------------------------------------------------------*/
char* concat(const char *s1, const char *s2);
char* concat2(const char *s1, const char *s2, const char *s3);
char* concat3(const char *s1, const char *s2, const char *s3, const char *s4);
char* concat4(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5);
char *concat_strings(int numberOfStr, char *basis, ...);

/* ------------------------------------------------------------*/
/* ARRAY INITIALIZATION                                        */
/* ------------------------------------------------------------*/
void initialize_array_char(int nbins, char *array, int value);
void initialize_array_int(int nbins, int *array, int value);
void initialize_array_int32_t(int nbins, int32_t *array, int value);
void initialize_array_long_int(int nbins, long int *array, int value);
void initialize_array_float(int nbins, float *array, float value);
void initialize_array_double(int nbins, double *array, double value);

void initialize_array_int_pointer(int nbins, int **array);
void initialize_array_float_pointer(int nbins, float **array);
void initialize_array_double_pointer(int nbins, double **array);

/* ------------------------------------------------------------*/
/* ARRAY ALLOCATION                                            */
/* ------------------------------------------------------------*/
char* allocate_array_char(int length, char *arrayname);
int* allocate_array_int(int length, char *arrayname);
int32_t* allocate_array_int32_t(int length, char *arrayname);
long int* allocate_array_long_int(int length, char *arrayname);
float* allocate_array_float(int length, char *arrayname);
double* allocate_array_double(int length, char *arrayname);

int** allocate_array_int_pointer(int length, char *arrayname);
float** allocate_array_float_pointer(int length, char *arrayname);
double** allocate_array_double_pointer(int length, char *arrayname);

/* ------------------------------------------------------------*/
/* FREE ALLOCATED MEMORY                                       */
/* ------------------------------------------------------------*/
void free_multiple(int numberOfThingsToFree, ...);

/* ------------------------------------------------------------*/
/* FILE EXISTENCE                                              */
/* ------------------------------------------------------------*/
int file_exist(char *filename);
int directory_exist(char *dirname);
void *get_directory(char *dirname);

/* ------------------------------------------------------------*/
/* CHRAS & STRINGS                                             */
/* ------------------------------------------------------------*/

char *replace_char(char* str, char find, char replace);

#endif
