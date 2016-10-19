#define L_STR 1000
#define bit_shift 8
#define N_CONF 20




int configurate(double CONFIG[N_CONF], char ITEM[N_CONF][L_STR], int already_read[N_CONF], char * err_msg, char * addCONF, char * prob);

void display_config(double * CONFIG, char * prob);

int initialize(int nInitValue,  realArray InitValue[nInitValue], int sizeIinitValue[nInitValue], char * err_msg, char addInitValue[nInitValue][L_STR], char * prob);


//void file_write_trouble(int m, int K, runList * runhist, char * sol_name);

//void write_column(int m, double data[], char * source, char * sol_name);



int DATA_OUTPUT(char * err_msg, int nDATA, char addDATA[nDATA][L_STR], int sizeDATA[nDATA], double * DATA[nDATA], int idx[nDATA], int flagDATA[nDATA], char * directory);
