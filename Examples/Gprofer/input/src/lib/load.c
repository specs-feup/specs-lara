/* Implementation file for lib/load */

#include <errno.h>
#include <inttypes.h>
#include "load.h"
#include "matisse.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tensor.h"
#include "tensor_struct.h"


/**
 */
tensor_d* load_matrix_variable_2_td(const char * varname, int dim_0, int dim_1, tensor_d** restrict out)
{

   new_array_d_2(dim_0, dim_1, out);
   load_raw_data_from_file_td(varname, (uint32_t) dim_0 * (uint32_t) dim_1, *out);
   
   return *out;
}

/**
 */
void load_raw_data_from_file_td(const char * varname, uint32_t num_elements, tensor_d* out)
{
// void load_raw_data_from_file (const char* varname, int num_elements, double* out)

    size_t element_size = sizeof(double);
    size_t read_total_elements = 0;
    
    char* relative_filename = malloc(sizeof("data/.dat") + strlen(varname));
    sprintf(relative_filename, "data/%s.dat", varname);
    char* filename = get_absolute_filename(relative_filename);
    free(relative_filename);
    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
       fprintf(stderr, "Could not open file %s.\n", filename);
       exit(1);
    }
    
    if (out->data == NULL) {
       fprintf(stderr, "Output data is null.\n");
       exit(1);
    }
    
    while (read_total_elements < num_elements) {
        size_t elements_to_read = num_elements - read_total_elements;
        size_t max_read = 1024 / element_size;
        if (elements_to_read > max_read) {
            elements_to_read = max_read;
        }
        size_t read_elements = fread(out->data + read_total_elements, element_size, elements_to_read, file);
        if (read_elements == 0) {
            // We should use %zu instead of %d for size_t printfs, but unfortunately some compilers
            // (such as MinGW) do not support it.
            fprintf(stderr, "Could not load file %s. Read %d of %d:\n%s\n", filename, (int) read_total_elements, (int) num_elements, strerror(errno));
            exit(1);
        }
        read_total_elements += read_elements;
    }
    
    free(filename);
    fclose(file);
}

/**
 */
char * get_absolute_filename(char * filename)
{
//char* get_absolute_filename(const char* filename) {

#ifdef _WIN32
   HMODULE module = GetModuleHandleA(NULL);
   char* file_path = NULL;
   size_t size = 32;
   do {
      size *= 2;
      file_path = realloc(file_path, size);
   } while (GetModuleFileNameA(module, file_path, size) >= size);
   
   // file_path has the path including the executable name now.
   
   size_t index = strlen(file_path) - 2;
   while (file_path[index] != '\\' && file_path[index] != '/') {
      if (--index <= 0) {
         fprintf(stderr, "Could not get executable path");
         exit(1);
      }
   }
   
   file_path[index + 1] = '\0';
#else
// disabled warning: #warning Fallback absolute path - returning relative path
   char* file_path = malloc(3);
   strcpy(file_path, "./");
#endif

   if (file_path == NULL) {
      fprintf(stderr, "Could not get executable folder\n");
      exit(1);
   }
   char* absolute_filename = malloc(strlen(file_path) + strlen(filename) + 1);
   if (absolute_filename == NULL) {
      fprintf(stderr, "Could not allocate memory for absolute filename\n");
      exit(1);
   }
   strcpy(absolute_filename, file_path);
   strcat(absolute_filename, filename);
   free(file_path);

   return absolute_filename;}
