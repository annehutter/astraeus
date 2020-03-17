#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <stdint.h>
#include <dirent.h>

/* ------------------------------------------------------------*/
/* RANDOM NUMBERS                                              */
/* ------------------------------------------------------------*/
double randdouble()
{
    return rand()/((double)(RAND_MAX)+1);
} 

/* ------------------------------------------------------------*/
/* STRING ADDITION                                             */
/* ------------------------------------------------------------*/
char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1) + strlen(s2) + 1);//+1 for the null-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    
    return result;
}

char* concat2(const char *s1, const char *s2, const char *s3)
{
    char *result = malloc(strlen(s1) + strlen(s2) + strlen(s3) + 1);//+1 for the null-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    
    return result;
}

char* concat3(const char *s1, const char *s2, const char *s3, const char *s4)
{
    char *result = malloc(strlen(s1) + strlen(s2) + strlen(s3) + strlen(s4) + 1);//+1 for the null-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    strcat(result, s4);
    
    return result;
}

char* concat4(const char *s1, const char *s2, const char *s3, const char *s4, const char *s5)
{
    char *result = malloc(strlen(s1) + strlen(s2) + strlen(s3) + strlen(s4) + strlen(s5) + 1);//+1 for the null-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    strcat(result, s4);
    strcat(result, s5);
    
    return result;
}

char *concat_strings(int numberOfStr, char *basis, ...)
{
    char *ptr;
    char *name;
    size_t length = strlen(basis);
    int count=1;

    va_list args_list;
    va_start(args_list, basis);

    for(int i=0; i<numberOfStr-1; i++)
    {
        ptr = va_arg(args_list, char*);
        length += strlen(ptr);
    }
    va_end(args_list);
    if((name = malloc(length + 1))) //+1 for the null-terminator   
    {
        va_start(args_list, basis);
        strcpy(name, basis);

        for(int i=0; i<numberOfStr-1; i++)
        {
            ptr = va_arg(args_list, char*);
            strcat(name, ptr);
            count++;
        }
        va_end(args_list);
        return name;
    }
    else
    {
        printf("Could not allocate memory for the name with basis :\n %s\n", basis);
        exit(EXIT_FAILURE);
    }
}

/* ------------------------------------------------------------*/
/* ARRAY INITIALIZATION                                        */
/* ------------------------------------------------------------*/
void initialize_array_char(int nbins, char *array, int value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

void initialize_array_int(int nbins, int *array, int value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

void initialize_array_int32_t(int nbins, int32_t *array, int value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

void initialize_array_long_int(int nbins, long int *array, int value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

void initialize_array_float(int nbins, float *array, float value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

void initialize_array_double(int nbins, double *array, double value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

/* ------------------------------------------------------------*/
/* ARRAY ALLOCATION                                            */
/* ------------------------------------------------------------*/

char* allocate_array_char(int length, char *arrayname)
{
    char *tmp = malloc(length * sizeof(char));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_char: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_char(length, tmp, '\0');
    return tmp;
}

int* allocate_array_int(int length, char *arrayname)
{
    int *tmp = malloc(length * sizeof(int));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_int: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_int(length, tmp, 0.);
    return tmp;
}

int32_t* allocate_array_int32_t(int length, char *arrayname)
{
    int *tmp = malloc(length * sizeof(int32_t));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_int: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_int32_t(length, tmp, 0.);
    return tmp;
}

long int* allocate_array_long_int(int length, char *arrayname)
{
    long int *tmp = malloc(length * sizeof(long int));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_int: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_long_int(length, tmp, 0.);
    return tmp;
}

float* allocate_array_float(int length, char *arrayname)
{
    float *tmp = malloc(length * sizeof(float));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_float: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_float(length, tmp, 0.);
    return tmp;
}

double* allocate_array_double(int length, char *arrayname)
{
    double *tmp = malloc(length * sizeof(double));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_double: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_double(length, tmp, 0.);
    return tmp;
}

/* ------------------------------------------------------------*/
/* FREE ALLOCATED MEMORY                                       */
/* ------------------------------------------------------------*/

void free_multiple(int numberOfThingsToFree, ...)
{
    void *ptr;
    va_list args_list;
    va_start(args_list, numberOfThingsToFree);
    for(int i=0; i<numberOfThingsToFree; i++)
    {
        ptr = va_arg(args_list, void *);
        free(ptr);
    }
    va_end(args_list);
}

/* ------------------------------------------------------------*/
/* FILE EXISTENCE                                              */
/* ------------------------------------------------------------*/

int file_exist(char *filename)
{
    FILE *file;
    if((file = fopen (filename, "rt"))){
        fclose(file);
        return 1;
    }else{
        return 0;
    }
}

int directory_exist(char *dirname)
{
        DIR* dir = opendir(dirname);
        if(dir)
        {
                closedir(dir);
                return 1;
        }
        else
        {
                return 0;
        }
}

void *get_directory(char *filename)
{
        char *token = NULL;
        char *directory = NULL;
        size_t length;

        token = strrchr(filename, '/');

        if( filename == NULL )
        {
                printf("There is no file:'%s'\n", filename); /* You decide here */
        }
        
        if(token != NULL)
        {
            length = strlen(token);
            directory = malloc(length);
            memcpy(directory, token+1, length);
        }

        return directory;
}

/* ------------------------------------------------------------*/
/* CHRAS & STRINGS                                             */
/* ------------------------------------------------------------*/

char *replace_char(char* str, char find, char replace)
{
    char *current_pos = strchr(str,find);
    while (current_pos){
        *current_pos = replace;
        current_pos = strchr(current_pos,find);
    }
    return str;
}