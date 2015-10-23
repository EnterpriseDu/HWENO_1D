#ifndef N_CONF
#define N_CONF 9
#endif /* N_CONF */
#ifndef N_OPT
#define N_OPT 8
#endif /* N_OPT */




typedef struct runNode{
  int *trouble0, *trouble1;
  char RcstrState[6];
  double RcstrErr[12];
  double time[2];
  struct runNode * next;
} runNode;

typedef struct{
  runNode * head, * current, * tail;
  int length;
} runList;


void init_runList(runList * runhist);

void insert_runList(runList * runhist);

void locate_runList(int p, runList * runhist);

void delete_runList(runList * runhist);




void mem_release();



/* This function examine whether a string
 * represents a real number.
 * Transform the string represents a
 * negtive number into a string represents
 * the opposite one and return its' sign.
 * It returns 0 if the string do not
 * represents a real number.
 * In the result given by this function, there will
 * be only one 'e' in the string, and the
 * only position for '-' is behind 'e', and
 * there can be only one dot in the string
 * and the only position for it is before 'e'.
 */
int format_string(char * str);


/* this function transform a string
 * consisting '1', '2', ..., and '.'
 * into the real number it represent
 */
double str2num(char * number);


/* this function counts how many lines are there
 * the initial data file.
 */
int data_pre_read_line(FILE * fp, int err_code);


/*
 * This function tells us how many columns
 * are there in the file.
 */
int data_pre_read_column(FILE * fp, int err_code);


/* This function reads the initial data file
 * to generate the initial data.
 * It returns 0 if successfully read the file,
 * while returns the index of the wrong entry.
 *
 * num = the expected number of entries in the file.
 */
int data_read(FILE * fp, double * U, int num);


/* this function counts how many lines are there
 * the configuration file.
 */
int string_pre_read_line(FILE * fp, int err_code);


/*
 * This function tells us how many columns
 * are there in the file.
 */
int string_pre_read_column(FILE * fp, int err_code);

int string_read(FILE * fp, char * str_conf[], int const n_conf, int const L);




int vec_read(char * add, int COLUMN, int err_code);


/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables* of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
int initialize(char * engine, char * prob, char * addrho, char * addu, char * addp, char * addx, int adp);


/* This function read the configuration data file,
 * and store the configuration data in the array
 * "CONFIG".
 * CONFIG[0] is the constant of perfect gas
 * CONFIG[1] is the CFL number
 * CONFIG[2] is the largest value can be seen as zero
 * CONFIG[3] is the first limiter of the slope
 * CONFIG[4] is the second limiter of the slope
 * CONFIG[5] is the first parameter of the monitor function
 * CONFIG[6] is the second parameter of the monitor function
 * CONFIG[7] is the modifier of the mesh redistribution
 * CONFIG[8] is the tolerance of the mesh redistribution
 */
int configurate(double * CONFIG, char * engine, char * prob, char * add);


/* OPT[0] is the maximal step to compute.
 * OPT[1] is the time to stop the computation
 * OPT[2] is the switch of whether keep the inter-data during the computation
 * OPT[3] is the switch of wether use an adaptive mesh
 */
int optionize(double * OPT, char * engine, char * prob, char * add);






/*
 *
 *
 *
 *
 */
void file_write_data
(int m, int start, int N, double * data[], char * source, char * sol_name);


/*
 *
 *
 *
 *
 */
void file_write_log
(int m, int n, int N, double scaling, double * CONFIG, double * OPT, runList * runhist, char * scheme, char * prob, char * sol_name);

void file_write_trouble(int m, int K, runList * runhist, char * sol_name);

void write_column(int m, double data[], char * source, char * sol_name);
