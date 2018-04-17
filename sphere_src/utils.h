#ifndef logio_h
#include "logio.h"
#endif

int string_copy(char **var, char **val);

char *query_user(const char query_str[]);

int lowstr(char s[]);

void char_free(char **cptr);

float frand(unsigned int iseed);

float uniform_rand(float num);

float gauss_rand(float num);

int two_dim_FHT(float data_array[], int xnum, int ynum, int forward);

int fast_hartley_transform(float data_array[], float work_array[], int num_el, int do_forward_transform);

int is_pow_of_two(int inum);

void fat_init(void);

FILE *check_lun(int lun, int want_read, int *is_text);

int file_open(char read_write_append, char is_text, char *filename);

int file_close(int lun);

int file_rewind(int lun);

int file_check_eof(int lun);

int file_error(int lun);

void file_cat(void);

int file_search(FILE *fptr, char string[], size_t slen, int leave_start);

int skip_field(FILE *fptr);

int input_array(FILE *fptr, float *array, size_t num_el);

#ifndef sphere_h
struct Descriptor;
#endif

int user_printf(FILE *fptr, char fmt[], int nargs, struct Descriptor *args[]);
int fmt_read(FILE *fptr, char fmt[], int nargs, struct Descriptor *args[]);

int ask_user(char prompt[]);
char *prompt_user(char *prompt, char *defstr);

int indexx(int npts, float arrin[], int **indx);

void get_increments(int axis[3], int ndim[3], int add[3]);

int  match(char *regexp, char *string, int *was_error);

int fourier_series(float x_data[], float y_data[], int npts, float period,
float asum[], float bsum[], int max_order);

int fourier_series_value(float xval, float *yval, int differential_order,
float period, float amp[], float phase[], float filter[], int max_order);

int match(char *regexp, char *string, int *was_error);
