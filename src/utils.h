#ifndef UTILS_H
#define UTILS_H

#ifndef MAXLENGTH
#define MAXLENGTH 1024
#endif

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

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
void initialize_array_int(int nbins, int *array, int value);
void initialize_array_int32_t(int nbins, int32_t *array, int value);
void initialize_array_long_int(int nbins, long int *array, int value);
void initialize_array_float(int nbins, float *array, float value);
void initialize_array_double(int nbins, double *array, double value);
void initialize_2darray_double(int nrows, int ncols, double **array, double value);

/* ------------------------------------------------------------*/
/* ARRAY ALLOCATION                                            */
/* ------------------------------------------------------------*/
int* allocate_array_int(int length, char *arrayname);
int32_t* allocate_array_int32_t(int length, char *arrayname);
long int* allocate_array_long_int(int length, char *arrayname);
float* allocate_array_float(int length, char *arrayname);
double* allocate_array_double(int length, char *arrayname);
double** allocate_2darray_double(int nrows, int ncols, char *arrayname);

/* ------------------------------------------------------------*/
/* FILE EXISTENCE                                              */
/* ------------------------------------------------------------*/
int file_exist(char *filename);
int directory_exist(char *dirname);
char *get_directory(char *dirname);

char *strdup(const char *s);
char *strndup(const char *s, size_t n);

/* ------------------------------------------------------------*/
/* COPYING FILES                                               */
/* ------------------------------------------------------------*/
void copy_file(char *fileToCopy, char *copiedTarget);

#endif
