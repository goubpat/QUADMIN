/*  _       _
 * | | _^_ | |  Statistics  Statistique
 * |_| >_< |_|  Canada      Canada
 *
 *
 * NAME
 *	quadminz.c
 *
 * DESCRIPTION
 *	Main program for the quadratic minimization program
 *
 * HISTORY
 *	Created  ?              ?*	Modified R. Puchyr      2003-03-10
 *		To reduce the possibility of error of the result series
 *		being indavertently multiplied by the frequency, legacy
 *		sections of the code have been removed in the get_ser()
 *		and upd_ser() functions that attempted to determine the
 *		power and the annunal rates (arates) attribute of the
 *		series from the database so that the resulting series 
 *		could be adjusted. The annual rates was a feature in the
 *		NASH data base, FAME's predessecor, which is not kept
 *		in FAME. Therefore such checks, and thus, the code is
 *		not necessary.
 *	Modified R. Puchyr      2003-04-23
 *		Off-by-one array out of bounds bug found and corrected
 *		in round.c function distribround(). This function was
 *		correctly identifying the ith element to adjust by one,
 *		so the results would sum to the bench, but it was
 *		adjusting the i-1th element.
 *	Modified Kevin Fredrickson Oct 24, 2017
 *		Exits after 100 attempts to get input from the procedure
 *		program. At this point we are assuming there is a 
 *		disconnect between the fame code and this c code.
 *		Also, limiting this code to run on 4 CPU's only. This
 *		code takes over a CPU. This allows other processes
 *		timely response.
 */




#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <errno.h>

#include <hli.h>


#include <time.h>

#define	YES	1
#define	NO	0
#define MAX_FAME_NAME 130   /* At least twice 64 because users can input: database_name'series_name as input  */
#define SHORT_BUF_SIZE 950

#define LANG_ENG	0
#define LANG_FRA	1


/* these defines are for NA, NC or ND */

#define	MISSNC	-999999.9999
#define	MISSND	-999998.9999
#define	MISSNA	-999997.9999



/*
 * Strcutures
 */

typedef char bool;

struct s_ser_info
{
	char base[65];
	int  benchfreq;
	int  freq;
	char from[7];
	char to[7];
	int  fiscallag;
};

struct s_algo
{
	bool linked;
	char linkto[7];
	bool round;
	int  decs;
	bool prop;
	bool first;
	bool update;
	char updatefrom[7];
	bool mean;
	bool stock;
	bool zero;
};

struct s_reports
{
	char file_name[BUFSIZ];
	bool display;
	bool arates;
	bool diff;
	bool fact;
	bool gr;
	int  lag;
};

struct s_series
{
	char benchid[65];
	char distributorid[65];
	char targetid[65];
};

struct s_options
{
	struct s_ser_info ser_info;
	struct s_algo     algo;
	struct s_reports  reports;
	struct s_series   series;
};



/*
 * Prototypes
 */

int init_base(struct s_options *opt);
void init_ser_info(struct s_ser_info *pnt);
void init_algo(struct s_algo *pnt, struct s_ser_info *ser_pnt);
void init_reports(struct s_reports *pnt);
void send_warning(struct s_options *opt, int setnum, int langnum, int messnum, char **parm, int nb_parm);
void open_output_file(char *file_name);
int read_fame_line(struct s_options *opt, char[]);
int get_fame_input(struct s_options *opt, int *still_job);
void end_fame(void);
int benchmark(struct s_options *opt, char bfrom[], char bto[]);
void upd_ser(struct s_options *opt, double *trget);
int write_ser(char *base, char *targetid, char *from, char *to, int freq, double *target);
void prnt_warnings(double *dist, double *trget, int nbdist, struct s_options *opt, char *bfrom);
void prnt_w_mess(struct s_options *opt, int num, char *mess1, char *mess2, int nbmess);
void roundser(double *trget, double *bench, int *tau, int *kappa, int nbbench, int nbdist, struct s_options *opt, char bto[]);
void print_reports(double *bench, double *dist, double *trget, int nbdist, int nbbench, struct s_options *opt, int *tau, int *kappa, double *af);
void prnt_replace(char **parm, int setnum, int langnum, int messnum, char *title, int nb_parm);
void cal_tau_kappa(int *tau, int *kappa, struct s_options *options, char bfrom[], char bto[]);
void add_date(char date[], int freq, int val);
int cal_nb_points(char from[], char to[], int freq, int freq2);
void ret_dates(struct s_options *pnt, char bfrom[], char bto[]);
void stock_start(struct s_options *pnt, char bfrom[]);
int divide(int per, int freq);
int get_ser(struct s_options *options, double **bench, double **dist, char bfrom[], char bto[]);
int read_series(char *base_name, int freq, char *from, char *to, double *out, char *ser_name);
void cal_fac(double *result, double *trget, double *dist, int nbdist, char prop);
void send_error(struct s_options *opt, char *short_buf);

extern void benchmod(double *x, double *b, double *cor, double *y, int *tau, int *kappa, double *w, int *prop, int *diff, int *index, int tt, int mm);
extern void print_default(double *dist, double *trget, char from[], int freq, int benchfreq, int nbpoints, int ndecs, int div, char stock, char *prnt);
extern void print_fisc(double *dist, double *trget, int *tau, int *kappa, int nbpoint, int nbbench, int ndecs, int freq, int benchfreq, char from[], int div, char stock, char *prnt);
extern void prnt_data(char start[], int nbpoints, int freq, int nbdecs, double *series, char arates, char printsum);
extern void read_print(FILE *out, int set_num, int lang_num, int mess_num, char **parm, int nb_parm, char *out_mess);

extern void firstdiff(int *nbdist, int *lag, double *dist, double *result);
extern void percent(int *nbdist, int *lag, double *dist, double *result);

extern char title[];
extern char *lookup_message(int setnum, int lang_num, int messnum);
extern char* replace(char *message, int nbtokens, char **tokens);

int distribround(double *in, double sum, int nvalue, int ndec, double *out);


int lang = 0;
int MAXTRY = 1000;
int workkey;
FILE *tables;
double mistt[3];



main()
{
	int  loop_ctr;
	char sys_cmd[36];	
	char pid[16];
	struct s_options options;
	char bfrom[7];
	char bto[7];
	int  still_job = 1;

	//	/**********
 //       * Restrict the CPU's on which this process can run.
	//**********/
	//strcpy(sys_cmd, "/bin/taskset -pc 4,5,6,7 ");
	///*snprintf(pid, 16, "%d", getpid());*/
	//strcat(sys_cmd, pid);
	//strcat(sys_cmd, " >/dev/null 2>&1");
	//system(sys_cmd);
	//*/

	/**********
	* initialisations
	**********/

	lang = LANG_ENG;

	init_base(&options);
	init_ser_info(&options.ser_info);
	init_algo(&options.algo, &options.ser_info);
	init_reports(&options.reports);

	/**********
	* The process is executed until the still job pointer is set to
	* false
	**********/
	loop_ctr = 0;
	while (still_job)
	{
		loop_ctr += 1;
		if (loop_ctr > 100) {
			/* At this point believe there is a disconnect between fame code and
			   this server code. Therefore we will just exit. */
			break;
		}
		if (!get_fame_input(&options, &still_job))
			continue;

		loop_ctr = 0;

		if (still_job == 0)
			break;

		if (options.reports.display)
			open_output_file(options.reports.file_name);

		ret_dates(&options, bfrom, bto);

		benchmark(&options, bfrom, bto);
	}
	end_fame();
}

/**********
 *
 * int init_base(struct s_options *opt)
 *
 * - Initialize the necessary chli function to interact with Fame
 * - Open workdatabase for process
 * - Translate missing values table
 *
 * returns: 1 if everything o.k.
 * Otherwise exit the program
 **********/

int init_base(struct s_options *opt)
{
	char short_buf[SHORT_BUF_SIZE];
	int status;

	cfmini(&status);

	if (status != HSUCC)
	{
		exit(-1);
	}

	cfmopwk(&status, &workkey);

	if (status != HSUCC)
	{
		if (lang == LANG_FRA)
			sprintf(short_buf, "Le Program ecrit en C n'a pu ouvrir la base de donnees WORK. Aucun traitement effectue");
		else
			sprintf(short_buf, "The C Program could not open the work database. No processing done");

		send_error(opt, short_buf);
		end_fame();
		exit(-1);
	}

	cfmspm(&status, MISSNC, MISSND, MISSNA, mistt);

	if (status != HSUCC)
	{
		if (lang == LANG_FRA)
			sprintf(short_buf, "Le Program ecrit en C n'a pu ouvrir obtenir les tables de translation");
		else
			sprintf(short_buf, "The C Program could not not set translation table.");

		send_error(opt, short_buf);
	}

	return(1);
}



/**********
 *
 * void init_ser_info(struct s_ser_info *pnt)
 *
 * set default values for series information
 *
 **********/

void init_ser_info(struct s_ser_info *pnt)
{
	strcpy(pnt->base, "");
	pnt->benchfreq = 1;
	pnt->freq = 4;
	pnt->from[0] = '\0';
	pnt->to[0] = '\0';
	pnt->fiscallag = 0;
}



/**********
 *
 * void init_algo(struct s_algo *pnt, struct s_ser_info *ser_pnt)
 *
 * set default values for algorithm information
 *
 **********/

void init_algo(struct s_algo *pnt, struct s_ser_info *ser_pnt)
{
	pnt->linked = NO;
	strcpy(pnt->linkto, ser_pnt->from);
	pnt->round  = NO;
	pnt->decs = 1;
	pnt->prop   = YES;
	pnt->first  = YES;
	pnt->update = NO;
	strcpy(pnt->updatefrom, ser_pnt->from);
	pnt->mean   = NO;
	pnt->stock  = NO;
}



/**********
 *
 * void init_reports(struct s_reports *pnt, struct s_ser_info *info);
 *
 * set default values for reports options
 *
 **********/

void init_reports(struct s_reports *pnt)
{
	strcpy(pnt->file_name, "");
	pnt->display = NO;
	pnt->arates = NO;
	pnt->diff = NO;
	pnt->fact = NO;
	pnt->gr = NO;
	pnt->lag = 1;
}



/**********
 *
 * int get_fame_input(struct s_options *opt, int *still_job)
 *
 * after calling the function which read input from the Fame procedure,
 * the function put the option into the data structure.
 *
 *
 *    return   1: everything o.k.
 *             0: else.
 *
 *********/


int get_fame_input(struct s_options *opt, int *still_job)
{
	static char input_line[SHORT_BUF_SIZE];
	int c = 0;

	do
	{
		if (c++ == MAXTRY)
		{
			*still_job = 0;
			if (lang == LANG_FRA)
				sprintf(input_line, "Le Program ecrit en C a rencontre une boucle infinie a la fonction get_fame_input");
			else
				sprintf(input_line, "The C program has encountered an infinite loop in the get_fame_input function ");

			send_error(opt,input_line);
			return(0);
			break;
		}

		if (!read_fame_line(opt,input_line))
		{
			return(0);
			break;
		}

		if (strncmp(input_line,"Q_BENCHFREQ",11) == 0)
		{
			opt->ser_info.benchfreq = atoi(&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_FREQ",6) == 0)
		{
			opt->ser_info.freq = atoi(&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_FROM",6) == 0)
		{
			strcpy(opt->ser_info.from,&input_line[20]);
			strcpy(opt->algo.linkto,&input_line[20]);
			strcpy(opt->algo.updatefrom,&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_TO",4) == 0)
		{
			strcpy(opt->ser_info.to,&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_BENCHID",9) == 0)
		{
			strcpy(opt->series.benchid,&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_DISTRIBUTORID",15) == 0)
		{
			strcpy(opt->series.distributorid,&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_TARGETID",10) == 0)
		{
			strcpy(opt->series.targetid,&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_BASE",6) == 0)
		{
			strcpy(opt->ser_info.base,&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_FISCALLAG",11)  == 0)
		{
			opt->ser_info.fiscallag = atoi(&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_LINKED",8) == 0)
		{
			opt->algo.linked = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_ROUND",7) == 0)
		{
			opt->algo.round = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_DECS",6) == 0)
		{
			opt->algo.decs = atoi(&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_PROP",6) == 0)
		{
			opt->algo.prop = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_FIRST",7) == 0)
		{
			opt->algo.first = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_UPDATE",8) == 0)
		{
			if (strncmp(input_line,"Q_UPDATEFROM",12) == 0)
			{
				strcpy(opt->algo.updatefrom,&input_line[20]);
				continue;
			}
			else
			{
				opt->algo.update = (input_line[20] == 'Y');
				continue;
			}
		}

		if (strncmp(input_line,"Q_MEAN",6) == 0)
		{
			opt->algo.mean = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_STOCK",7) == 0)
		{
			opt->algo.stock = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_ZERO",6) == 0)
		{
			opt->algo.zero = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_DISPLAY",9) == 0)
		{
			opt->reports.display = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"QUAD_TABLE_OUTPUT",17) == 0)
		{
			strcpy(opt->reports.file_name, &input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_ARATES",8) == 0)
		{
			opt->reports.arates = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_DIFF",6) == 0)
		{
			opt->reports.diff = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_FACT",6) == 0)
		{
			opt->reports.fact = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_GR",4) == 0)
		{
			opt->reports.gr = (input_line[20] == 'Y');
			continue;
		}

		if (strncmp(input_line,"Q_LAG",5) == 0)
		{
			opt->reports.lag = atoi(&input_line[20]);
			continue;
		}

		if (strncmp(input_line,"Q_LANG",6) == 0)
		{
			lang = (strcmp(&input_line[20],"FRA") == 0 ? LANG_FRA : lang);
			continue;
		}

		if (strncmp(input_line,"END",3) == 0)
		{
			*still_job = 0;
			break;
		}


	} while ((strncmp(input_line,"Q_ARG_PAST",10) != 0) && (strncmp(input_line,"END",3) != 0));

	return(1);
}



/**********
 *
 * read_fame_line(struct s_options *opt,char input_line[])
 *
 * Get input sent by the fame procedure
 *
 *
 *    return   1: everything o.k.
 *             0: else.
 *
 *********/

int read_fame_line(struct s_options *opt,char input_line[])
{
	char short_buf[SHORT_BUF_SIZE];
	int status;

	cfmsinp(&status, input_line);

	if (status != HSUCC)
	{
		if (lang == LANG_FRA)
			sprintf(short_buf, "Le Program ecrit en C n'a pu lire tous les intrants provenant de Fame.");
		else
			sprintf(short_buf, "The C program could not read some input from Fame.");

		send_error(opt, short_buf);
		return(0);
	}

	return(1);
}



/**********
 *
 * void open_output_file(char *file_name)
 *
 * opens file_name for append.
 *
 * If the file does not contain a dot, the extension ".txt" will be added to
 * the file name.
 *
 **********/

void open_output_file(char *file_name)
{
	char fname[BUFSIZ];

	if (strcmp(file_name, "") != 0)
	{
		strcpy(fname, file_name);

		if (strchr(fname, '.') == NULL)
			strcat(fname, ".txt");

		if ((tables = fopen(fname, "a")) == NULL)
			tables = stdout;
	}
	else
		tables = stdout;

}



/**********
 *
 * ret_dates(struct s_options *pnt, char bfrom[], char bto[])
 *
 * This routine calculates the retrieval dates for benchmark and distributor.
 *
 * From dates
 *
 * For benchmarks, it compares with the from date and calculates the closest
 * benchmark start date greater or equal to the from date.
 *
 * To dates
 *
 * Calculates the greatest benchmark end date possible that still covers
 * the distributor data.  It's done considering the fiscal adjustment.
 *
 * bfrom and bto are the start and end dates for benchmarks.
 *
 **********/

void ret_dates(struct s_options *pnt, char bfrom[], char bto[])
{
	int  per, year;

	/**********
	* from dates
	**********/

	strcpy(bfrom, pnt->ser_info.from);
	strcpy(bto, pnt->ser_info.to);
	add_date(bfrom, pnt->ser_info.freq, -pnt->ser_info.fiscallag);
	per = atoi(bfrom+4);

	/**********
	* If the series are stock, the way to calculate retrieval is
	* different.
	**********/

	if (pnt->algo.stock)
	{
		stock_start(pnt, bfrom);
	}
	else
	{
		while (!divide(per,pnt->ser_info.benchfreq))
		{
			add_date(bfrom, pnt->ser_info.freq,1);
			per = atoi(bfrom+4);
		}
		if (pnt->ser_info.benchfreq == 4)
		{
			per = (per + 2) / 3;
			sprintf(bfrom+4, "%02.2d", per);
		}
	}

	/**********
	* be sure not to update earlier than the linkto point.
	**********/

	if (strcmp(pnt->algo.updatefrom, pnt->algo.linkto) <= 0 && pnt->algo.linked)
	{
		 strcpy(pnt->algo.updatefrom, pnt->algo.linkto);
		 add_date(pnt->algo.updatefrom, pnt->ser_info.freq, 1);
	}

	/**********
	* to date
	**********/

	add_date(bto, pnt->ser_info.freq, -(pnt->ser_info.fiscallag));
	per = atoi(bto+4);
	year = atoi(bto) - per;
	if (pnt->ser_info.benchfreq == 1)
	{
		if (per == pnt->ser_info.freq)
			year = year + 1;         /* y per 1 */
		else
			year = year - 100 + 1;   /* y - 1 per 1 */
	}
	else
	{
		 while (!(per == 3 || per == 6 || per == 9 || per == 12))
		 {
			 add_date(bto, pnt->ser_info.freq, (int)-1);
			 per = atoi(bto + 4);
			 year = atoi(bto) - per;
		 }

		 year += per/3;
	}

	sprintf(bto, "%.6d", year);
}


/**********
 *
 * void add_date(char date[], int freq, int val)
 *
 * This procedure performs addition or substractions on dates.
 * The specified date is a string and is added up with the int val.
 * The date has to be in the valid format of yyyypp.  Where yyyy is a year
 * and pp is a period that can range from one to frequency.  It adds or
 * substract but always keeps the pp in the range from 1 to frequency (
 * frequency beeing 04 for quaterly and 12 for monthly).
 *
 * note: This procedure is not meant to add big number for val.  If
 * so, it would be much more efficient to rewrite this procudure using
 * modulos.
 *
 **********/

void add_date(char date[], int freq, int val)
{
	int per;
	int year;

	per = atoi(date+4);
	year = atoi(date) - per;
	per = per + val;

	while (per > freq)
	{
		per -= freq;
		year += 100;
	}

	while (per < 1)
	{
		per += freq;
		year -= 100;
	}

	year += per;
	sprintf(date, "%.6d", year);
}



/**********
 *
 * void stock_start(struct s_options* pnt, char bfrom[]);
 *
 * calculates the start date for retrieval of benchmark when the
 * series is stock.  When coming in, bfrom = from.
 * the date calculated is :
 * when benchfreq == 1 ->  sameyear + first period.
 * when benchfreq == 4 ->  sameyear + if from period =  1  2  3 -> 1
 *                                                      4  5  6 -> 2
 *                                                      7  8  9 -> 3
 *                                                     10 11 12 -> 4
 **********/

void stock_start(struct s_options *pnt, char bfrom[])
{
	int per;

	per = atoi(bfrom+4);
	if (pnt->ser_info.benchfreq == 1)
		per = 1;
	else
		per = ((per - 1) / 3) + 1;

	sprintf(bfrom+4, "%02.2d", per);
}



/**********
 *
 * int divide(int per, int freq)
 *
 * checks to see if a number in a given period is compatible with
 * the frequency.  This has to be done when we copy
 * a date from one frequency to another.
 * The possible frequencies are 1 and 4 (benchmark frequency)
 * if freq = 1   per can range from 1 to 12 (freq = 4 or 12 for dist)
 * if freq = 4   per can range from 1 to 12 (freq = 12 for dist)
 *
 **********/

int divide(int per, int freq)
{
	if (freq == 1)
	{
		if (per == 1)
			return(1);
		else
			return(0);
	}
	else
	{
		if (per == 1 || per == 4 || per == 7 || per == 10)
			return(1);
		else
			return(0);
	}
}



/**********
 *
 * int  benchmark(struct s_options *opt, char bfrom[], char bto[])
 *
 * - Calls the function to read the series
 * - Allocates all the space needed for the calculations
 * - Calls the function to calculate the reference points
 * - Calls the function to execute the benchmarking algorithm
 * - if needed, calls the function to round the series
 *
 *    return   1: everything o.k.
 *             0: else.
 *
 **********/

int benchmark(struct s_options *opt, char bfrom[], char bto[])
{
	double *bench;
	double *dist;
	double *cor;
	double *trget;
	double *weights;
	int *tau;
	int *kappa;
	int i, nbdist, nbbench, j;
	int prop, diff, index;
	char short_buf[SHORT_BUF_SIZE];

	prop = (opt->algo.prop  ? 0 : 1);
	diff = (opt->algo.first ? 1 : 2);
	index = (opt->algo.mean  ? 1 : 0);


	/**********
	* read the series
	**********/

	if (!get_ser(opt, &bench, &dist, bfrom, bto))
		return(0);

	nbdist = cal_nb_points(opt->ser_info.from, opt->ser_info.to, opt->ser_info.freq, opt->ser_info.freq);
	nbbench = cal_nb_points(bfrom, bto, opt->ser_info.benchfreq, opt->ser_info.benchfreq);
	if (opt->algo.linked)
		nbbench++;

	tau     =    (int *)malloc(nbbench * sizeof(int));
	kappa   =    (int *)malloc(nbbench * sizeof(int));
	cor     = (double *)malloc(nbdist  * sizeof(double));
	trget   = (double *)malloc(nbdist  * sizeof(double));
	weights = (double *)malloc((nbdist+1)  * sizeof(double));


	if (!(tau && kappa && cor && trget && weights))
	{
		if (lang == LANG_FRA)
			sprintf(short_buf, "Le Program ecrit en C n'a pu allouer assez de memoire. Essayer des series plus courtes");
		else
			sprintf(short_buf, "The C program could not allocate memory. You might want to try smaller series");

		send_error(opt, short_buf);
		return(0);
	}

	/**********
	* calculate reference points
	**********/

	cal_tau_kappa(tau, kappa, opt, bfrom, bto);

	for (i = 0; i < (nbdist+1); i++)
		weights[i] = 1.0;

	/**********
	*
	* There is one error to check in this procedure
	* it occurs when  tau[0] = tau[1]  &
	*                 stock = y
	*
	**********/

	if (tau[0] == tau[1] && opt->algo.stock)
	{
		prnt_w_mess(opt, 0, "", "", 0);  /* The index number is 10 , incremented in prnt_w_mess  */

		for (i = 2; i < nbbench; i++)
		{
			bench[i-1] = bench[i];
			tau[i-1]   = tau[i];
			kappa[i-1] = kappa[i];
		}

		nbbench--;
	}

	/**********
	* benchmarking algorithm
	**********/

	(void)benchmod(dist, trget, cor, bench, tau, kappa, weights, &prop, &diff, &index, nbdist, nbbench);

	/**********
	* round if needed
	*********/

	if (opt->algo.round)
	{
		roundser(trget, bench, tau, kappa, nbbench, nbdist, opt, bto);
	}

	/**********
	* If option zero, set related target values to zero if benchmark is zero
	**********/

	if (opt->algo.zero)
	{
		for (i = 0; i < nbbench ; i++)
		{
			if (bench[i] == 0)
			{
				for (j = tau[i] -1 ; j <= kappa[i] -1 ; j++)
					trget[j] = 0;
			}
		}
	}


	/**********
	* print warning messages
	**********/

	prnt_warnings(dist, trget, nbdist, opt, bfrom);

	/**********
	* update if needed
	**********/

	if (opt->algo.update)
		upd_ser(opt, trget);

	/**********
	* print the reports if needed
	**********/

	if (opt->reports.display)
		print_reports(bench, dist, trget, nbdist, nbbench, opt, tau, kappa, cor);

	free(bench);
	free(dist);
	free(tau);
	free(kappa);
	free(cor);
	free(trget);
	free(weights);

	return(1);
}



/**********
 *
 * int  get_ser(struct s_options *options, double **bench, double **dist,
 *              char bfrom[], char bto[])
 *
 * Function that reads the series from the work database.
 *
 * returns: 1 if everything o.k.
 *          0 else
 *
 **********/

int get_ser(struct s_options *options, double **bench, double **dist, char bfrom[], char bto[])
{
	struct s_ser_info *pnt;
	int nbbench, nbdist;
	int start;
	char base[MAX_FAME_NAME];
	char freq[3];
	char minimum[7];
	char cont;
	char bench_bool;
	char less;
	char benchid[MAX_FAME_NAME];
	char trgetid[MAX_FAME_NAME];
	char distid[MAX_FAME_NAME];
	char short_buf[SHORT_BUF_SIZE];

	strcpy(benchid, options->series.benchid);
	strcpy(trgetid, options->series.targetid);
	strcpy(distid, options->series.distributorid);


	start = 0;
	pnt = &options->ser_info;
	strcpy(base, pnt->base);
	strcpy(freq, "* ");


	/**********
	* calculate number of point to retrieve for benchmarks and dist.
	**********/

	nbbench = cal_nb_points(bfrom, bto, pnt->benchfreq, pnt->benchfreq);
	nbdist = cal_nb_points(pnt->from, pnt->to, pnt->freq, pnt->freq);

	/**********
	* if linked, need to retrieve one more benchmark.
	**********/

	if (options->algo.linked)
	{
		nbbench++;
		start = 1;
	}

	/**********
	* Allocate space.
	**********/

	*bench = (double *) malloc((size_t)(nbbench * sizeof(double)));
	*dist  = (double *) malloc((size_t)(nbdist  * sizeof(double)));


	if (!(*bench && *dist))
	{
		if (lang == LANG_FRA)
			sprintf(short_buf, "Le Program ecrit en C n'a pu allouer assez de memoire. Essayer des series plus courtes");
		else
			sprintf(short_buf, "The C program could not allocate enough memory. You might want to try smaller series");

		send_error(options, short_buf);
		return(0);
	}

	/**********
	* If linked, first benchmark datapoint is retrieved from target
	* series.  The date to retrieve from is in linkto.
	**********/

	if (options->algo.linked)
	{
		sprintf(freq,"%02.2d", pnt->freq);
		if (read_series(base, pnt->freq, options->algo.linkto, options->algo.linkto, *bench, trgetid) != 1)
		{
			if (lang == LANG_FRA)
				sprintf(short_buf, "Le Program ecrit en C n'a pu lire la serie cible");
			else
				sprintf(short_buf, "The C Program could not read the target series");

			send_error(options, short_buf);
			return(0);
		}
	}

	/**********
	* get benchmark data
	**********/

	cont = (char)2;
	strcpy(minimum, bfrom);
	add_date(minimum, pnt->benchfreq, 1);
	sprintf(freq,"%02.2d", pnt->benchfreq);
	bench_bool = (char)1;
	less = (char)0;
	while (cont == 2)
	{
		if (strcmp(minimum, bto) > 0)
			bench_bool = (char)0;

		cont = read_series(base, pnt->benchfreq, bfrom, bto, &((*bench)[start]), benchid);

		if (cont == 2)
		{
			add_date(bto, pnt->benchfreq, -1);
			less = 1;
		}

		if (strcmp(bfrom, bto) >= 0 || cont == 0)
		{
			if (lang == LANG_FRA)
				sprintf(short_buf, "Le Program ecrit en C n'a pas pu lire la serie jalon");
			else
				sprintf(short_buf, "The C Program could not read the benchmark series");

			send_error(options, short_buf);
			return(0);
		}
	}

	if (less)
		prnt_w_mess(options, 7, "", "", 0);	 /*  Missing data end of bechmark.  The index number is 17 , incremented in prnt_w_mess  */

	if (!bench_bool)
		prnt_w_mess(options, 11, "", "", 0);	 /* Only on benchmark used.  The index number is 21 , incremented in prnt_w_mess  */

	/**********
	* get distributor data
	**********/

	sprintf(freq,"%02.2d", pnt->freq);
	if (read_series(base, pnt->freq, pnt->from, pnt->to, *dist, distid) != 1)
	{
		if (lang == LANG_FRA)
			sprintf(short_buf, "Le Program ecrit en C n'a pas pu lire la serie distributrice");
		else
			sprintf(short_buf, "The C Program could not read the distributor series");

		send_error(options, short_buf);
		return(0);
	}

	return(1);
}



/**********
 *
 * int cal_nb_points(char from[], char to[], int freq, int freq2)
 *
 * calculate the number of point in between 2 dates
 * the 2 extremities are included
 *
 **********/

int cal_nb_points(char from[], char to[], int freq, int freq2)
{
	int syear, sdate, eyear, edate;
	int nbpoints;

	sdate = atoi(from+4);
	edate = atoi(to+4);
	syear = (atoi(from) - sdate)/100;
	eyear = (atoi(to) - edate)/100;

	if (freq == 4 && freq2 == 12)
		sdate = (sdate * 3) - 2;
	if (freq == 12 && freq2 == 4)
		edate = (edate * 3) - 2;

	nbpoints = (eyear - syear) * freq + edate - sdate + 1;
	return(nbpoints);
}



/**********
 *
 * int read_series(char *base_name, int freq, char *from, char *to,
 *                 double *out, char *ser_name)
 *
 * to read a series
 *
 * the series is read through the use of cfmfame function. It is done by
 * passing commands to fame and producing a temporary series (tmp) in the
 * work database.  The work database is already opened (init_base function).
 *
 * returns:	1 	if everything o.k.
 *		0 	else.
 *		2	Missing values found
 *
 **********/

int read_series(char *base_name, int freq, char *from, char *to, double *out, char *ser_name)
{
	int ret_val;
	int numobs;
	int syear, sprd;
	int eyear, eprd;
	int status;
	int tfreq;
	int range[3];
	int i;
	bool missing;
	char fame_cmd[BUFSIZ];
	char tmp_ser_name[MAX_FAME_NAME];

	ret_val = 1;
	switch (freq)
	{
		case 1:
			tfreq = HANDEC;
			break;

		case 4:
			tfreq = HQTDEC;
			break;

		case 12:
			tfreq = HMONTH;
			break;
	}

	numobs = -1;
	sprd  = atoi(from+4);
	eprd = atoi(to+4);
	syear = atoi(from);
	eyear = atoi(to);
	syear = (syear - sprd) / 100;
	eyear = (eyear - eprd) / 100;


	sprintf(fame_cmd, "date %04.4d to %04.04d; copy <overwrite on> %s as Q_TMP_SER to WORK", syear, eyear, ser_name);

	cfmfame(&status, fame_cmd);

	if (status != HSUCC)
		return(0);


	cfmsrng(&status, tfreq, &syear, &sprd, &eyear, &eprd, range, &numobs);

	if (status != HSUCC)
		return(0);


	strcpy(tmp_ser_name, "Q_TMP_SER");

	cfmrrng(&status, workkey, tmp_ser_name, range, out, HTMIS, mistt);


	if (status != HSUCC)
	{
		return(0);
	}
	else
	{
		missing = 0;
		for (i = 0; i < numobs && missing == 0; i++)
			missing = (out[i] == MISSNC || out[i] == MISSND || out[i] == MISSNA);

		if (missing)
			ret_val = 2;
	}

	return(ret_val);
}



/**********
 *
 * void cal_tau_kappa(int *tau, int *kappa, struct s_options *options,
 *                    char bfrom[], char bto[])
 *
 * calculate the tau and kappa arrays needed for the benchmod routine.
 * tau and kappa represent the starting and ending of the reference
 * period of the benchmarks over the series to adjust.  This procedure
 * takes care of fiscal year adjustment, linked options, and benchmarks
 * that are not available for the whole span of the distributor.
 *
 **********/

void cal_tau_kappa(int *tau, int *kappa, struct s_options *options, char bfrom[], char bto[])
{
	int i, j, inc;
	int start;
	int nbpoints;
	int gap, gap2;

	/**********
	* if the linked option is specified, put the linkto point at
	* the start of the tau and kappa tables.  Linkto is always equal to
	* the from date except for one case:  When fiscallag is < 0 then
	* the from date can be recalculated but the linkto stays the same.
	**********/

	start = 0;
	if (options->algo.linked)
	{
		start = 1;
		if (strcmp(options->algo.linkto, options->ser_info.from) == 0)
		{
			tau[0] = 1;    /* linkto = from date => first point */
			kappa[0] = 1;
		}
		else   /* fiscallag < 0 */
		{
			gap2 = cal_nb_points(options->ser_info.from, options->algo.linkto, options->ser_info.freq, options->ser_info.freq);
			tau[0] = gap2;
			kappa[0] = gap2;
		}
	}

	/**********
	* calculate number of benchmark points
	**********/

	nbpoints = cal_nb_points(bfrom, bto, options->ser_info.benchfreq, options->ser_info.benchfreq);

	/**********
	* find out how many points there is between start of
	* distributor and start of benchmarks
	**********/

	gap = cal_nb_points(options->ser_info.from, bfrom, options->ser_info.freq, options->ser_info.benchfreq) - 1;
	if (options->ser_info.benchfreq == 4)  /* freq = 12 */
		inc = 3;
	else
		inc = (options->ser_info.freq == 12 ? 12 : 4);

	/**********
	* calculate first reference (not linked) taking care of difference
	* in start points between dist. and benchmark (gap)
	**********/

	tau[start] = 1 + gap;
	kappa[start] = inc + gap;

	/**********
	* if stock option tau = kappa
	**********/
	if (options->algo.stock)
		tau[start] = kappa[start];

	/**********
	* recalculate first point if fiscallag != 0. shift over.
	**********/

	if (options->ser_info.fiscallag)
	{
		tau[start]   += options->ser_info.fiscallag;
		kappa[start] += options->ser_info.fiscallag;
	}

	/**********
	* calculate rest of tables
	* note that if linked = y the size of the array is
	* nbpoints + 1 so we don't do bounds violations.
	**********/

	for (i = start + 1, j = 1; j < nbpoints; j++, i++)
	{
		tau[i]   = tau[i-1]   + inc;
		kappa[i] = kappa[i-1] + inc;
	}
}




/**********
 *
 * void roundser(double *trget, double *bench, int *tau, int *kappa,
 *               int nbbench, int nbdist, struct s_options *opt)
 *
 * This procedure will round the target series up to a number of
 * decimal points.  It does so by calling the distribround function
 * with the appropriate figures.
 * We loop and call distribround with every benchmark and corresponding
 * target values.  For target values that are not covered by benchmark
 * we pass the Yearly (or quaterly) sum.
 *
 * tau and kappa give the reference periods of the benchmark on the
 * targets. tau = start period. kappa = end period.
 *
 * Modified:
 * Daniel Sheehy  Dec 1997
 * Adeed rounding to incomplete periods and not covered by benchmark.
 *
 **********/

void roundser(double *trget, double *bench, int *tau, int *kappa, int nbbench, int nbdist, struct s_options *opt, char bto[])
{
	int start;
	int initial_start;
	int nvalue, ndecs;
	int i,j;
	int len;
	int remain;
	unsigned int stcrounding;
	double result[12];
	double sum;

	start = 0;
	if (opt->algo.linked)
		start = 1;
	initial_start = start;
	stcrounding = 0;
	ndecs = opt->algo.decs;



	/**********
	* periods that are covered by benchmarks
	**********/

	for (i = start; i < nbbench; i++)
	{
		start = tau[i] -1;
		nvalue = kappa[i] - start;

		distribround(&trget[start], bench[i], nvalue, ndecs, result);
		for (j = 0; j < nvalue; j++)
		{
			trget[start+j] = result[j];
		}
	}

	/**********
	* periods not complete before benchmarks
	**********/

	remain = (tau[initial_start] -1) % nvalue;
	if (remain != 0)
	{
		for (j = 0, sum = 0; j < remain; j++)
			sum += (double)trget[initial_start+j];

		distribround(&trget[initial_start], sum, nvalue, ndecs, result);
		for (j = 0; j < nvalue; j++)
			trget[initial_start+j] = result[j];
	}

	/**********
	* periods that are not covered by benchmarks
	*********/

	len = (nbdist - (start+j)) / nvalue;
	remain = (nbdist - (start+j)) % nvalue;
	if (remain != 0)
		len++;

	start = start + nvalue;

	for (i = 0 ; i < len ; i++, start +=  nvalue)
	{
		if (i == len -1  && remain != 0)
			nvalue = remain;

		for (j = 0, sum = 0; j < nvalue; j++)
			sum += (double)trget[start+j];

		distribround(&trget[start], sum, nvalue, ndecs, result);
		for (j = 0; j < nvalue; j++)
			trget[start+j] = result[j];
	}
}



/**********
 *
 * print_reports(double *bench, double *dist, double *trget, int nbdist,
 *               int nbbench, struct s_options *opt, int *tau, int *kappa,
 *               double *af)
 *
 * Procedure to print the reports.  One report is  printed by
 * default and there 3 other kind of reports that can be printed if
 * asked for.
 *
 **********/

void print_reports(double *bench, double *dist, double *trget, int nbdist, int nbbench, struct s_options *opt, int *tau, int *kappa, double *af)
{
	double *result;
	int lag;
	int benchfreq, freq;
	int ndec;
	int i, div, start;
	int setnum;
	char from[7];
	char title[BUFSIZ];
	char **parameter;
	char prnt;


	/**********
	* this first part is to calculate which number to pass to the
	* report function.  This number will divide the sum of the targets.
	* the sum have to be divided for printing purposes only if
	* mean = y.
	**********/

	div = 1;
	if (opt->algo.mean)
	{
		if (opt->ser_info.benchfreq == 1)
			div = (opt->ser_info.freq == 12 ? 12 : 4);
		else
			div = 3;
	}

	parameter = (char **)malloc(1 * sizeof(char *));
	setnum = 4;

	start = 0;
	strcpy(from, opt->ser_info.from);
	ndec = opt->algo.decs;
	freq = opt->ser_info.freq;
	benchfreq = opt->ser_info.benchfreq;
	lag = 1;
	result = (double *)malloc(nbdist * sizeof(double));

	/**********
	* print default report
	**********/

	fprintf(tables, " \n");

	prnt_replace(parameter, setnum, lang, 7, title, 0);
	fprintf(tables, " BENCHID       = %s\n", opt->series.benchid);
	fprintf(tables, " DISTRIBUTORID = %s\n", opt->series.distributorid);
	fprintf(tables, " TARGETID      = %s\n\n", opt->series.targetid);
	prnt_replace(parameter, setnum, lang, 16, title, 0);
	print_default(dist, trget, from, freq, benchfreq, nbdist, ndec, div,
				  opt->algo.stock, &prnt);
	if (prnt)
		prnt_replace(parameter, setnum, lang, 13, title, 0);

	fprintf(tables, "\n\n\n");

	if (opt->ser_info.fiscallag)
	{
		if (opt->algo.linked)
		{
			/**********
			* if linked, first benchmark is only one point which
			* is the linkto point.
			* - don't send this benchmark to printing module (start = 1)
			* - number of benchmark is one less (nbbench--)
			**********/

			start = 1;
			nbbench--;
		}

		prnt_replace(parameter, setnum, lang, 8, title, 0);
		print_fisc(dist, trget, &tau[start], &kappa[start], nbdist, nbbench, ndec, freq, benchfreq, from, div, opt->algo.stock, &prnt);

		if (prnt)
			prnt_replace(parameter, setnum, lang, 14, title, 0);

		fprintf(tables, "\n\n\n");

		if (opt->algo.linked)
			nbbench++;
	}

	parameter[0] = opt->series.distributorid;
	prnt_replace(parameter, setnum, lang, 6, title, 1);

	prnt_data(from, nbdist, freq, ndec, dist, opt->reports.arates, (char)1);
	fprintf(tables, "\n\n\n");

	parameter[0] = opt->series.targetid;
	prnt_replace(parameter, setnum, lang, 9, title, 1);

	prnt_data(from, nbdist, freq, ndec, trget, opt->reports.arates, (char)1);

	/**********
	* print adjustment factors if asked.
	**********/

	if (opt->reports.fact)
	{
		fprintf(tables, "\n\n\n");
		cal_fac(af, trget, dist, nbdist, opt->algo.prop);
		for (i = 0; i < nbdist; i++)
			result[i] = af[i];

		prnt_replace(parameter, setnum, lang, 10, title, 0);

		prnt_data(from, nbdist, freq, ndec, result, (char)NO, (char)NO);
	}

	/**********
	* print first differences if asked
	**********/

	if (opt->reports.diff)
	{
		fprintf(tables, "\n\n\n");
		prnt_replace(parameter, setnum, lang, 11, title, 0);
		parameter[0] = opt->series.distributorid;
		prnt_replace(parameter, setnum, lang, 18, title, 1);

		firstdiff(&nbdist, &lag, dist, result);
		prnt_data(from, nbdist, freq, ndec, result, (char)NO, (char)NO);

		fprintf(tables, "\n\n\n");
		prnt_replace(parameter, setnum, lang, 11, title, 0);

		parameter[0] = opt->series.targetid;
		prnt_replace(parameter, setnum, lang, 19, title, 1);

		firstdiff(&nbdist, &lag, trget, result);
		prnt_data(from, nbdist, freq, ndec, result, (char)NO, (char)NO);
	}

	/**********
	* print growth rate if asked
	**********/

	if (opt->reports.gr)
	{
		lag = opt->reports.lag;

		fprintf(tables, "\n\n\n");
		prnt_replace(parameter, setnum, lang, 12, title, 0);

		parameter[0] = opt->series.distributorid;
		prnt_replace(parameter, setnum, lang, 18, title, 1);

		percent(&nbdist, &lag, dist, result);
		prnt_data(from, nbdist, freq, 2, result, (char)NO, (char)NO);

		fprintf(tables, "\n\n\n");
		prnt_replace(parameter, setnum, lang, 12, title, 0);

		parameter[0] = opt->series.targetid;
		prnt_replace(parameter, setnum, lang, 19, title, 1);

		percent(&nbdist, &lag, trget, result);
		prnt_data(from, nbdist, freq, 2, result, (char)NO, (char)NO);
		printf("\n");
	}

	free(parameter);
	free(result);
	fflush(tables);
}



/**********
 *
 * void prnt_replace(char **parm, int setnum, int langnum, int messnum, char *title,
 *                   int nb_parm)
 *
 *
 **********/

void prnt_replace(char **parm, int setnum, int langnum, int messnum, char *title, int nb_parm)
{
	read_print(tables, setnum, langnum, messnum, parm, nb_parm, title);
}



/**********
 *
 * void cal_fac(double *result, double *trget, double *dist, int nbdist, char prop);
 *
 *   Calculate adjustment factors
 *   if prop = y then divide
 *   else substract;
 *
 **********/

void cal_fac(double *result, double *trget, double *dist, int nbdist, char prop)
{
	int i;

	if (prop)
		for (i = 0; i < nbdist; i++)
			result[i] = trget[i] / dist[i];
	else
		for (i = 0; i < nbdist; i++)
			result[i] = trget[i] - dist[i];
}



/**********
 *
 * void upd_ser(struct s_options *opt, double *trget)
 *
 * procedure to update the target series.
 * we update the series from the updatefrom date.
 *
 **********/

void upd_ser(struct s_options *opt, double *trget)
{
	int start, end;
	char base[MAX_FAME_NAME];
	char trgetid[MAX_FAME_NAME];
	char short_buf[SHORT_BUF_SIZE];

	/**********
	 *
	 * calculate points from which to start updating and
	 * how many points to update.
	 *
	 */
	start = cal_nb_points(opt->ser_info.from, opt->algo.updatefrom, opt->ser_info.freq, opt->ser_info.freq) - 1;
	end   = cal_nb_points(opt->ser_info.from, opt->ser_info.to, opt->ser_info.freq, opt->ser_info.freq);

	strcpy(base, opt->ser_info.base);
	strcpy(trgetid, opt->series.targetid);


	/**********
	 *
	 * Write the series data
	 *
	 */
	if (!write_ser(base, trgetid, opt->algo.updatefrom, opt->ser_info.to, opt->ser_info.freq, trget+start))
	{
		if (lang == LANG_FRA)
			sprintf(short_buf, "Le Program ecrit en C n'a pu mettre a jour la serie cible");
		else
			sprintf(short_buf, "The C Program could not update the target series");

		send_error(opt, short_buf);
	}
}



/**********
 *
 * int write_ser(char *base, char *targetid, char *from, char *to, double *target)
 *
 * to write data into the work database (The temporary target series).
 * The real target series will be updated in the Fame procedure.
 *
 * returns:	1 if everything o.k.
 *		    0 else
 *
 **********/

int write_ser(char *base, char *targetid, char *from, char *to, int freq, double *target)
{
	int numobs;
	int syear, sprd;
	int eyear, eprd;
	int status;
	int tfreq;
	int range[3];
	char tmp_series_name[] = "Q_TMP_UPDATED_SER";
	int i;

	switch (freq)
	{
		case 4:
			tfreq = HQTDEC;
			break;

		case 12:
			tfreq = HMONTH;
			break;
	}

	numobs = -1;
	sprd  = atoi(from+4);
	eprd = atoi(to+4);
	syear = atoi(from);
	eyear = atoi(to);
	syear = (syear - sprd) / 100;
	eyear = (eyear - eprd) / 100;

	cfmsrng(&status, tfreq, &syear, &sprd, &eyear, &eprd, range, &numobs);
	if (status != HSUCC)
	{
		return(0);
	}

	cfmwrng(&status, workkey, tmp_series_name, range, target, HNTMIS, mistt);
	if (status != HSUCC)
	{
		return(0);
	}

	return(1);
}



/**********
 *
 * void prnt_warnings(double *dist, double *trget, int nbdist,
 *                    struct s_options *opt, char *bfrom);
 *
 * To check for a few problems cases and call functions to print
 *  some warning messages if they occur.
 *
 **********/

void prnt_warnings(double *dist, double *trget, int nbdist, struct s_options *opt, char *bfrom)
{
	int i;
	int per;
	bool isneg;
	char optimal[7];
	char upd_optimal[7];

	strcpy(optimal, bfrom);
	if (opt->ser_info.benchfreq == 4)
	{
		per = atoi(optimal+4);
		per = per * 3;
		per = per - 2;
		sprintf(optimal+4, "%02.2d", per);
	}

	add_date(optimal, opt->ser_info.freq, opt->ser_info.fiscallag);
	strcpy(upd_optimal, optimal);
	add_date(optimal, opt->ser_info.freq, -1);

	/**********
	* check for unsatisfied benchmarks and movement discontinuity
	**********/
	if (opt->algo.linked)
	{
		if (opt->algo.stock)
		{ 							/*calculate optimal position*/
		}

		if (strcmp(optimal, opt->algo.linkto))
		{
			prnt_w_mess(opt, 1, optimal, "", 1);
			if (!opt->algo.stock)
				prnt_w_mess(opt, 2, opt->algo.updatefrom, optimal, 2);
		}
	}

	if (strcmp(upd_optimal, opt->algo.updatefrom) && opt->algo.update)
	{
		prnt_w_mess(opt, 3, upd_optimal, "", 1);
		if (!opt->algo.stock)
			prnt_w_mess(opt, 4, opt->algo.updatefrom, upd_optimal, 2);
	}

	/**********
	* check target for negative values
	**********/
	isneg = NO;
	for (i = 0; i < nbdist && !isneg; i++)
		isneg = (trget[i] < 0.0);

	if (isneg)
		prnt_w_mess(opt, 5, "", "", 0);

	/**********
	* check distributor for values <= 0
	**********/
	isneg = NO;

	if (opt->algo.zero)
	{
		for (i = 0; i < nbdist && !isneg; i++)
			isneg = (dist[i] < 0.0);
	}
	else
 	{
		for (i = 0; i < nbdist && !isneg; i++)
			isneg = (dist[i] <= 0.0);
	}

	if (isneg)
		prnt_w_mess(opt, 6, "", "", 0);

}



/**********
 *
 * void prnt_w_mess(struct s_options *opt, int num, char *mess1, char *mess2, int nbmess);
 *
 * Evaluates parameters and call send_warning function.
 *
 **********/

void prnt_w_mess(struct s_options *opt, int num, char *mess1, char *mess2, int nbmess)
{
	int setnum;
	int messnum;
	char **parameter;

	parameter = (char **)malloc(2 * sizeof(char *));
	parameter[0] = mess1;
	parameter[1] = mess2;
	setnum = 1;

	messnum = num + 10;    /* num = 1 -> messnum = 11 */
	send_warning(opt, setnum, lang, messnum, parameter, nbmess);
	free(parameter);
}



/**********
 *
 * void send_warning(struct s_options *opt, int setnum, int langnum,
 *                 int messnum, char **parm, int nb_parm);
 *
 * Send a warning message through Fame warning channel.
 *
 *********/

void send_warning(struct s_options *opt, int setnum, int langnum, int messnum, char **parm, int nb_parm)
{
	int status;
	char short_buf[SHORT_BUF_SIZE];
	char fame_cmd[BUFSIZ];

	strcpy(short_buf, lookup_message(setnum, langnum, messnum));

	if (nb_parm)
		replace(short_buf, nb_parm, parm);

	if (lang == LANG_FRA)
		sprintf(fame_cmd, "signal warning : \"MESSAGE QUADMIN pour les series %s, %s, %s : \"", opt->series.benchid,opt->series.distributorid,opt->series.targetid);
	else
		sprintf(fame_cmd, "signal warning : \"QUADMIN MESSAGE for %s, %s, %s: \"", opt->series.benchid,opt->series.distributorid,opt->series.targetid);

	cfmfame(&status,fame_cmd);

	sprintf(fame_cmd, "signal warning : \"%s\"",short_buf);
	cfmfame(&status,fame_cmd);
}



/**********
 *
 * void send_error(struct s_options *opt, char *short_buf);
 *
 * Send an error message through Fame error channel.
 *
 *********/

void send_error(struct s_options *opt, char *short_buf)
{
	int status;
	char fame_cmd[BUFSIZ];

	if (lang == LANG_FRA)
		sprintf(fame_cmd, "signal continue : \"MESSAGE QUADMIN pour les series %s, %s, %s: \" +newline + \"%s\" +newline", opt->series.benchid,opt->series.distributorid,opt->series.targetid,short_buf);
	else
		sprintf(fame_cmd, "signal continue : \"QUADMIN MESSAGE for %s, %s, %s: \" +newline + \"%s\" +newline", opt->series.benchid,opt->series.distributorid,opt->series.targetid,short_buf);

	cfmfame(&status,fame_cmd);
}



/**********
 *
 * void end_fame(void)
 *
 * Terminate interaction with Fame
 *
 **********/

void end_fame(void)
{
	int status;

	cfmfin(&status);
}



void benchmod(double *x, double *b, double *cor, double *y,
	int *tau, int *kappa, double *w, int *prop,
	int *diff, int *index, int tt, int mm);

void build_qinvw(double *x2, double *x, double *rquinv, int *tau,
	int tt, int *kappa, double *w, int prop,
	double *qinvw, int mm);

void build_wqinvw(int *tau, int *kappa, int mm, double *qinvw,
	double *wqinvw, double *w, int tt);

void cal_discrep(double mm, int *tau, int *kappa,
	double *add_disc, double *pro_disc, double *y,
	double *x, double *w, int index);

double *cal_inv(int dim, double *mat);
double *cal_inv2(int dim, double *mat);
void matmult(double *aa, double *bb, double *cc, int rowb, int colc, int colb);
void apply_corr(int tt, double *b, double *x, double *cor, int prop);
void modif_corr(int *kappa, double *cor, int tt, double *b, double *x, int mm, int prop);
double sumit(double *x, int nbelem);
void send_out_of_mem(void);



/* int nbweights = 120; */

/**********
 *
 * this program does the quadratic minimization algorithm
 * it is a c implementation of a sas program given to us by
 * Pierre Cholette. 
 *
 * This version is also fully compatible between micro and hp.
 * there is only one line to add or remove depending when working
 * on hp or micro respectively.  It is the first line of all wich
 * is :
 *    #pragma prototypes on
 *
 *
 **********/

/*********
 *
 *  main routine that calls all the others to do the quadmin
 *  algorithm
 *
 *********/

void benchmod(double *x, double *b, double *cor, double *y,
	int *tau, int *kappa, double *w, int *prop,
	int *diff, int *index, int tt, int mm)
{
	double  *qinvw;
	double  *wqinvw;
	double  *wqinvw2;
	double  *rquinv;
	double  *x2;
	double  *add_disc;
	double  *pro_disc;
	double  *invy;
	size_t   size;

	if (*prop != 1)
		*prop = 0;
	if (*diff != 2)
		*diff = 1;
	if (*index != 1)                  /* index = 1 for index series   */
		*index = 0;

	size     = (size_t)sizeof(double);
	qinvw    = (double *)malloc(size * (size_t)(tt * mm));
	wqinvw   = (double *)malloc(size * (size_t)(mm * mm));
	rquinv   = (double *)malloc(size * (size_t)(tt));
	x2       = (double *)malloc(size * (size_t)(tt));
	add_disc = (double *)malloc(size * (size_t)(mm));
	pro_disc = (double *)malloc(size * (size_t)(mm));
	invy     = (double *)malloc(size * (size_t)(mm));

	if (!(qinvw && wqinvw && rquinv && cor && x2 && add_disc && pro_disc && invy))
		send_out_of_mem();

	build_qinvw(x2, x, rquinv, tau, tt, kappa, w, *prop, qinvw, mm);

	build_wqinvw(tau, kappa, mm, qinvw, wqinvw, w, tt);

	cal_discrep(mm, tau, kappa, add_disc, pro_disc, y, x, w, *index);

	wqinvw2 = cal_inv2(mm, wqinvw);

	matmult(invy, wqinvw2, add_disc, mm, (int)1, mm);

	matmult(cor, qinvw, invy, tt, (int)1, mm);

	apply_corr(tt, b, x, cor, *prop);
	if (*diff == 2)
		modif_corr(kappa, cor, tt, b, x, mm, *prop);

	free(qinvw);
	free(wqinvw);
/**********
 * next line is needed only if the matrix inversion routine is
 * cal_inv.  If the other one is picked out (cal_inv2),
 * comment out the free statement.
 **********/
/*  free(wqinvw2);   */
	free(rquinv);
	free(x2);
	free(add_disc);
	free(pro_disc);
	free(invy);
}

/*********
 *
 * construction of the matrix qinvw.
 *
 * this function takes by far the most time to run of all the
 * functions in this program.  Because of this, it has been rearranged
 * to take a few shortcuts but this is at the expense of readability.
 * Mainly, what was done was to take as much calculation as possible out of
 * the loops.  In order to do so, it was necessary to use some temporary
 * variables which are all declared in the second part of the declaration block.
 *
 **********/

void build_qinvw(double *x2, double *x, double *rquinv, int *tau,
	int tt, int *kappa, double *w, int prop,
	double *qinvw, int mm)
{
	int r,c,m,k;
	int expo;
	double tw1;
	double t1;
	double xbar;
	double nperm;
	double rho;

	double temp;
	double tpow;
	double tdiv;
	double *trquinv;
	double *tw;
	double *tx2;

	xbar = sumit(x, tt) / tt;
	rho  = 0.99999999;

	for (r = 0; r < tt; r++)
		x2[r] = (prop == 0) ? x[r] : 1;

	for (r = 0; r < tt; r++)
	{
		trquinv = rquinv; 
		tx2 = x2;
		tdiv = x2[r] / xbar;
		for (c = 0; c < tt; c++, tx2++, trquinv++)
		{
			expo = abs(c-r);
			tpow = pow(rho, (double)expo);
			*trquinv = tpow * tdiv * *tx2;
		}

		tw1 = 0;
		for (m = 0; m < mm; m++)
		{
			t1 = tau[m] - 1;
			nperm = kappa[m] - t1;
			temp = 0;
			trquinv = &rquinv[(int)t1];
			tw = &w[(int)tw1];
			for (k = 0; k < nperm; k++, trquinv++, tw++)
				temp += *trquinv * *tw;
			tw1 = tw1 + nperm;
			*qinvw = temp;
			qinvw++;
		}
	}
}

/**********
 *
 * construction of the matrix wqinvw
 *
 **********/

void build_wqinvw(int *tau, int *kappa, int mm, double *qinvw, double *wqinvw, double *w, int tt)
{
	double tw1;
	double t1;
	double nperr;
	double temp;
	int    r,c,k;

	tw1 = 0;
	for (r = 0; r < mm; r++)
	{
		t1 = tau[r] - 1;
		nperr = kappa[r] - t1;

		for (c = 0; c < mm; c++)
		{
			temp = 0;
			for (k = 0; k < nperr; k++)
				temp += w[(int)(tw1+k)] * qinvw[(int)(((t1+k)*mm)+c)];

			*wqinvw = temp;
			wqinvw++;
		}

		tw1 += nperr;
	}
}

/**********
 *
 * calculates the annual discrepancies
 *
 **********/

void cal_discrep(double mm, int *tau, int *kappa,
	double *add_disc, double *pro_disc, double *y,
	double *x, double *w, int index)
{
	int r,k;
	double t1;
	double tw1;
	double temp;
	double nperr;
	double yr;
	double sumw;
	double *tw;
	double *tx;

	tw1 = 0;
	for (r = 0; r < mm; r++)
	{
		t1 = tau[r] - 1;
		nperr = kappa[r] - t1;
		temp = 0;
		sumw = 0;
		tw = &w[(int)tw1];
		tx = &x[(int)t1];

		for (k = 0; k < nperr; k++, tw++, tx++)
		{
			sumw += *tw;
			temp += *tw * *tx;
		}

		yr = (index == 1) ? y[r] * sumw : y[r];
		add_disc[r] = yr - temp;
		pro_disc[r] = yr / temp;
		tw1 += nperr;
	}
}

/**********
 *
 * calculates the inverse of a matrix. 
 *
 * algorithm:
 *     we proceed by the 'triangularisation' of the matrix mat.
 *     we do the same operations on the identity matrix (inv) then on
 *     the mat matrix.  when the 'triangularisation' of mat is done
 *     the inv matrix as become the inverse of the matrix mat.
 *
 * Note: Not used. cal_inv2 is used instead. 
 * 
 **********/

double *cal_inv(int dim, double *mat)
{
	int row,col,row2;
	int     add;
	double *inv;
	double *tempmat;
	double  pivot;

	inv = (double *)malloc((size_t)(dim * dim * sizeof(double)));
	tempmat = (double *)malloc((size_t)(dim * dim * sizeof(double)));

	if (! (inv && tempmat))
		send_out_of_mem();

	for (row = 0; row < dim; row++)
	{
		add = row*dim;
		for (col = 0; col < dim; col++)
		{
			tempmat[add+col] = mat[add+col];
			inv[add+col] = (row == col) ? 1 : 0;
		}
	}

	for (row = 0; row < dim; row++)
	{
		add = row*dim;
		pivot = tempmat[add+row];
		for (col = 0; col < dim; col++)
		{
			tempmat[add+col] /= pivot;
			inv[add+col] /= pivot;
		}
		for (row2 = 0; row2 < dim; row2++)
		{
			if (row2 != row)
			{
				pivot = -tempmat[row2*dim+row];
				for (col = 0; col < dim; col++)
				{
					tempmat[row2*dim+col] += pivot * tempmat[add+col];
					inv[row2*dim+col] += pivot * inv[add+col];
				}
			}
		}
	}

	free(tempmat);
	return(inv);
}

/**********
 * 
 * calculates the inverse of a matrix.
 * 
 * Finds the invrese of a positive definite symmetrix matrix by the
 * pivotal condensation method. this routine was furnished by
 * P. Cholette it does the same as the function above but is supposed
 * to be more precise.
 *
 **********/

double *cal_inv2(int dim, double *mat)
{
	int     icol, l, l1;
	int     add,add2;
	int     mult;
	double  pivot;
	double  t;

	for (icol = 0; icol < dim; icol++)
	{
		mult = 1;
		add = icol*dim;
		pivot = mat[add+icol];
		mat[add+icol] = 1;
		if (pivot < 0)
			mult = -1;
		if (fabs(pivot) < 1.0e-20)
			pivot = 1.0e-20 * mult;
		pivot = 1 / pivot;

		for (l = 0; l < dim; l++)
			mat[add+l] *=pivot;

		for (l1 = 0; l1 < dim; l1++)
		{
			if (l1 != icol)
			{
				add2 = l1*dim;
				t = mat[add2+icol];
				mat[add2+icol] = 0;
				for (l = 0; l < dim; l++)
					mat[add2+l] = mat[add2+l] - mat[add+l] * t;
			}
		}
	}

	return(mat);
}

/*********
 *
 * this routine does a matrix multiplication of matrix bb and cc and
 * puts the result of the multiplication into matrix aa.
 *
 * the matrixes are passed in parameters as pointers to doubles so
 * we have to calculate each address of the matrixes.
 *
 * matrix bb as dimensions rowb * colb
 * matrix cc as dimensions colb * colc   (rowc as to be = to colb)
 * matrix aa as dimensions rowb * colc
 *
 **********/

void matmult(double *aa, double *bb, double *cc, int rowb, int colc, int colb)
{
	int r,c,k;
	int add;
	double sum;

	for (r = 0; r < rowb; r++)
	{
		add = r*colb;
		for (c = 0; c < colc; c++)
		{
			sum = 0.0;
			for (k = 0; k < colb; k++)
				sum += bb[add+k] * cc[k*colc+c];

			*aa = sum;
			aa++;
		}
	}
}

/**********
 *
 * applies the corrections to obtain the benchmarked series
 *
 **********/

void apply_corr(int tt, double *b, double *x, double *cor, int prop)
{
	int r;

	for (r = 0; r < tt; r++)
	{
		b[r] = x[r] + cor[r];
		if (prop == 0)
			cor[r] = b[r] / x[r];
	}
}

/**********
 *
 * modifies the correction at the end, if second differences chosen.
 *
 **********/

void modif_corr(int *kappa, double *cor, int tt, double *b, double *x, int mm, int prop)
{
	int t2;
	double last_dif;
	int   r;

	t2 = kappa[mm-1];
	if (t2 < tt)
	{
		last_dif = cor[(int)t2-1] - cor[(int)(t2-2)];
		for (r = t2; r < tt; r++)
		{
			cor[r] = cor[r-1] + last_dif;
			if (prop == 0)
				b[r] = x[r] * cor[r];
			else
				b[r] = x[r] + cor[r];
		}
	}
}


  
/**********
 *
 * sums the elements of the array x.
 *
 **********/

double sumit(double *x, int nbelem)
{
	int i;
	double result;

	for (result = 0, i = 0; i < nbelem; i++)
		result += x[i];
	return(result);
}


/**********
 * 
 * void send_out_of_mem(void);
 *
 * Send error message nor enough memory
 *
 **********/

void send_out_of_mem(void)
{
	int  status;
	char fame_cmd[BUFSIZ];

	sprintf(fame_cmd, "signal continue : \"QUADMIN ERROR: The C program could not allocate enough memory for calculation \" +newline"); 
	cfmfame(&status, fame_cmd);
	exit(-1);
} 



/*
C      PURPOSE:  CALCULATES THE PERCENTAGE CHANGE OF A SERIES
C       ON THE BASIS OF EITHER PERIOD TO PERIOD OR PERIOD OVER THE
C       SAME PERIOD OF THE PRECEEDING YEAR, DEPENDING UPON THE
C       UNDERLYING SEASONALITY OF THE TIME SERIES DATA.  THE VALUE
C       0.00005 HAS BEEN ARBITRARILY CHOSEN TO DENOTE DIVISION BY ZERO.
*/

void percent(int *nopts, int *lag, double *data, double *rslt)
{
	int j, l;
	l = *lag;

	for (j=0; j < *nopts; j++)
	{
		rslt[j]=0.0;
	}

	if (*lag < *nopts)
	{
		for (j=l; j< *nopts; j++)
		{
			if (fabs(data[j - *lag]) > 0.00005)
			{
				rslt[j]=((data[j]-data[j - *lag])/data[j - *lag])*100;
			}
		}
	}
}



/*
C       PURPOSE:  CALCULATES FIRST DIFFERENCES EITHER PERIOD TO PERIOD
C       OR PERIOD OVER THE SAME PERIOD OF THE PRECEEDING YEAR, DEPENDING
C       UPON THE SEASONALITY OF THE UNDERLYING TIME SERIES.
C
C
C       IN THE EVENT THAT THE NUMBER OF OBSERVATIONS IS LESS
C       THAN THE LAG THEN THE RESULT VECTOR WILL BE FILLED
C       WITH ZERO VALUES.
*/

void firstdiff(int *nopts, int *lag, double *data, double *rslt)
{
	int l, i;

	l=*lag;
	for (i=0; i<*nopts; i++)
	{
		rslt[i]=0;
	}

	for (i=l; i< *nopts; i++)
	{
		rslt[i] = data[i] - data[i-*lag];
	}
}


#define PARM      "\"!\""

void read_print(FILE *out, int setnum, int langnum, int messnum, char **parm, int nb_parm, char * buf);
char *replace(char *message, int nbtokens, char **tokens);

extern char *lookup_message(int setnum, int lang_num, int messnum);
extern int errno;



/**********
 *
 * void read_print(FILE *out, int setnum, int langnum,
 *                 int messnum, char **parm, int nb_parm, char * buf);
 *
 * 1- read a message from the message file
 * 2- if needed, replace parameters in the message
 * 3- print the resulting string to the specied output file
 *
 **********/

void read_print(FILE *out, int setnum, int langnum, int messnum, char **parm, int nb_parm, char * buf)
{
	strcpy(buf, lookup_message(setnum, langnum, messnum));

	if (nb_parm > 0)
		replace(buf, nb_parm, parm);

	fprintf(out, "%s\n", buf); /* newline added by R Puchyr 2003-03-12 */
}



/**********
 *
 * char *replace(char *message, int nbtokens, char **tokens)
 *
 * takes a messages and replaces each occurence of the string PARM by
 * a token.  Be sure there is enough space to put everything in the
 * string message.
 * returns a pointer to the changed string.
 *
 **********/

char *replace(char *message, int nbtokens, char **tokens)
{
	int   nbfound;
	char *rest;
	char  rest2[BUFSIZ];

	nbfound = 0;

	while ( (rest = strstr(message, PARM)) && (nbfound < nbtokens) )
	{
		rest[0] = '\0';
		rest += strlen(PARM);
		strcpy(rest2, rest);
		strcat(message, tokens[nbfound]);
		strcat(message, rest2);
		nbfound++;
	}

	return(message);
}




char *lookup_message(int setnum, int langnum, int messnum)
{
	static char ret[1024];

	switch (setnum)
	{
		case 1:
			switch (langnum)
			{
				case 1:
					switch (messnum)
					{
						case  1:
							strcpy(ret, "ERREUR 1: ARGUMENT INVALIDE \"!\" ASSOCIE AU MOT CLEF \"!\".");
							break;

						case  2:
							strcpy(ret, "ERREUR 2: CONFLIT ENTRE BENCHFREQ ET FREQ.");
							break;

						case  3:
							strcpy(ret, "ERREUR 3: DUREE SPECIFIEE EST TROP COURTE.");
							break;

						case  4:
							strcpy(ret, "ERREUR 4: CONFLIT DE DATES ENTRE UPDATEFROM ET FROM OU TO.");
							break;

						case  5:
							strcpy(ret, "ERREUR 5: DISTRIBUTORID ET TARGETID NE PEUVENT AVOIR LE MEME NOM. VEUILLEZ RE-ENTRER LA VALEUR ASSOCIEE AU MOT CLEF TARGETID.");
							break;

						case  6:
							strcpy(ret, "ERREUR 6: COMMANDE INCONNUE \"!\" = \"!\".");
							break;

						case  7:
							strcpy(ret, "ERREUR 7: L'ARGUMENT SPECIFIE POUR FISCALLAG DOIT ETRE PLUS PETIT QUE FREQ.");
							break;

						case  8:
							strcpy(ret, "ERREUR 8: LA COMBINAISON D'OPTIONS: MEAN=Y ET STOCK=Y EST INVALIDE.");
							break;

						case  9:
							strcpy(ret, "ERREUR 9: VALEUR MANQUANTE POUR LE MOT-CLEF: \"!\".");
							break;

						case 10:
							strcpy(ret, "ERROR 10: LE MEME DISTRIBUTEUR NE PEUT ETRE EGAL A DEUX JALONS DIFFERENTS.");
							break;

						case 11:
							strcpy(ret, "IL EST POSSIBLE QU'IL Y AIT DISCONTINUITE DE MOUVEMENT ENTRE LA SERIE CALCULEE ET CELLE QUI LA PRECEDE SUR LA BASE. DATE DE DEBUT (FROM) RECOMMANDEE: \"!\".");
							break;

						case 12:
							strcpy(ret, "IL EST POSSIBLE QUE LE JALON COUVRANT LA PERIODE \"!\" NE SOIT PAS SATISFAIT. DATE DE DEBUT (FROM) RECOMMANDEE: \"!\".");
							break;

						case 13:
							strcpy(ret, "IL EST POSSIBLE QU'IL Y AIT DISCONTINUITE DE MOUVEMENT ENTRE LA SERIE MISE A JOUR ET LES CHIFFRES QUI PRECEDENT SUR LA BASE. DATE DE MISE A JOUR (UPDATEFROM) RECOMMANDEE: \"!\".");
							break;

						case 14:
							strcpy(ret, "IL EST POSSIBLE QUE LE JALON COUVRANT LA PERIODE \"!\" NE SOIT PAS SATISFAIT. DATE DE MISE A JOUR (UPDATEFROM) RECOMMANDEE: \"!\".");
							break;

						case 15:
							strcpy(ret, "LA SERIE CIBLE A DES VALEURS NEGATIVES.");
							break;

						case 16:
							strcpy(ret, "LA SERIE DISTRIBUTRICE A DES VALEURS PLUS PETITES OU EGALES A ZERO.");
							break;

						case 17:
							strcpy(ret, "VALEURS MANQUANTES A LA FIN DE LA SERIE JALON. DISTRIBUTEUR EXTRAPOLE.");
							break;

						case 18:
							strcpy(ret, "ERREUR 10:   LA SERIE CIBLE NE PEUT ETRE LUE");
							break;

						case 19:
							strcpy(ret, "ERREUR 11:   LA SERIE DISTRIBUTRICE NE PEUT ETRE LUE VEUILLEZ VERIFIER LES DATE AINSI QUE LA FREQUENCE ");
							break;

						case 20:
							strcpy(ret, "ERROR 10:   LA SERIE JALON NE PEUT ETRE LUE");
							break;

						case 21:
							strcpy(ret, "SEULEMENT UNE VALEUR JALON UTILISEE.");
							break;

						default:
							sprintf(ret, "No message for number s%dm%d", setnum, messnum);
							break;

					}

					break;

				case 2:
					switch (messnum)
					{
						case  1:
							strcpy(ret, "ÎØÈÁÊÀ 1: ÍÅÄÎÏÓÑÒÈÌÎÅ ÇÍÀ×ÅÍÈÅ \"!\" ÄËß ÊËÞ×ÅÂÎÃÎ ÑËÎÂÀ: \"!\".");
							break;

						case  2:
							strcpy(ret, "ERROR 2: ÍÅÑÎÎÒÂÅÒÑÒÂÈÅ ÌÅÆÄÓ BENCHFREQ È FREQ.");
							break;

						case  3:
							strcpy(ret, "ERROR 3: ÇÀÄÀÍÍÛÉ ÂÐÅÌÅÍÍÎÉ ÄÈÀÏÀÇÎÍ ÍÅÄÎÑÒÀÒÎ×ÍÎ ÂÅËÈÊ.");
							break;

						case  4:
							strcpy(ret, "ERROR 4: ÍÅÑÎÎÒÂÅÒÑÒÂÈÅ ÌÅÆÄÓ UPDATEFROM È ÄÀÒÀÌÈ FROM ÈËÈ TO.");
							break;

						case  5:
							strcpy(ret, "ERROR 5: DISTRIBUTORID È TARGETID ÍÅ ÌÎÃÓÒ ÈÌÅÒÜ ÎÄÍÎ È ÒÎ ÆÅ ÈÌß.  ÏÎÆÀËÓÉÑÒÀ, ÏÎÂÒÎÐÈÒÅ ÂÂÎÄ TARGETID.");
							break;

						case  6:
							strcpy(ret, "ERROR 6: ÍÅÏÀÑÏÎÇÍÀÍÍÀß ÏÀÐÀ KEYWORD=VALUE: \"!\" = \"!\".");
							break;

						case  7:
							strcpy(ret, "ERROR 7: ÇÀÄÀÍÍÎÅ ÇÍÀ×ÅÍÈÅ FISCALLAG ÄÎËÆÍÎ ÁÛÒÜ ÌÅÍÜØÅ FREQ.");
							break;

						case  8:
							strcpy(ret, "ERROR 8: ÍÅ ÄÎÏÓÑÊÀÅÒÑß ÎÄÍÎÂÐÅÌÅÍÍÎÅ ÇÀÄÀÍÈÅ MEAN=Y AND STOCK=Y.");
							break;

						case  9:
							strcpy(ret, "ERROR 9: ÎÒÑÓÒÑÒÂÓÞÙÅÅ ÇÍÀ×ÅÍÈÅ ÄËß ÊËÞ×ÅÂÎÃÎ ÑËÎÂÀ: \"!\".");
							break;

						case 10:
							strcpy(ret, "ÍÅÂÎÇÌÎÆÍÎ ÏÐÈÐÀÂÍßÒÜ ÎÄÈÍ È ÒÎÒ ÆÅ ÒÅÊÓÙÈÉ ÏÎÊÀÇÀÒÅËÜ ÄÂÓÌ ÐÀÇËÈ×ÍÛÌ ÊÎÍÒÐÎËÜÍÛÌ ÇÍÀ×ÅÍÈßÌ. ÏÅÐÂÎÅ ÈÇ ÍÈÕ ÈÑÏÎËÜÇÎÂÀÍÎ ÄËß ÓÂßÇÊÈ, ÂÒÎÐÎÅ ÏÐÎÈÃÍÎÐÈÐÎÂÀÍÎ.");
							break;

						case 11:
							strcpy(ret, "ÅÑÒÜ ÂÅÐÎßÒÍÎÑÒÜ ÐÀÇÐÛÂÀ Â ÄÈÍÀÌÈÊÅ ÌÅÆÄÓ ÊÎÐÐÅÊÒÈÐÓÅÌÛÌ Â ÄÀÍÍÛÉ ÌÎÌÅÍÒ ÐßÄÎÌ È ÏÐÅÄÛÄÓÙÈÌ ÐßÄÎÌ Â ÁÀÇÅ. ÐÅÊÎÌÅÍÄÓÅÒÑß ÍÀ×ÀÒÜ Ñ ÄÀÒÛ: \"!\".");
							break;

						case 12:
							strcpy(ret, "ÊÎÍÒÐÎËÜÍÛÉ ÏÎÊÀÇÀÒÅËÜ, ÎÕÂÀÒÛÂÀÞÙÈÉ ÏÅÐÈÎÄ \"!\" , ÂÅÐÎßÒÍÎ, ÍÅ ÑÎÃËÀÑÓÅÒÑß Ñ ÊÎÐÐÅÊÒÈÐÓÅÌÛÌ Â ÄÀÍÍÛÉ ÌÎÌÅÍÒ ÐßÄÎÌ. ÐÅÊÎÌÅÍÄÓÅÒÑß ÍÀ×ÀÒÜ Ñ ÄÀÒÛ: \"!\".");
							break;

						case 13:
							strcpy(ret, "ÅÑÒÜ ÂÅÐÎßÒÍÎÑÒÜ ÐÀÇÐÛÂÀ Â ÄÈÍÀÌÈÊÅ ÌÅÆÄÓ ÊÎÐÐÅÊÒÈÐÓÅÌÛÌ Â ÄÀÍÍÛÉ ÌÎÌÅÍÒ ÐßÄÎÌ È ÏÐÅÄÛÄÓÙÈÌ ÐßÄÎÌ Â ÁÀÇÅ. ÐÅÊÎÌÅÍÄÓÅÒÑß UPDATEFROM ÍÀ×ÈÍÀß Ñ ÄÀÒÛ: \"!\".");
							break;

						case 14:
							strcpy(ret, "ÊÎÍÒÐÎËÜÍÛÉ ÏÎÊÀÇÀÒÅËÜ, ÎÕÂÀÒÛÂÀÞÙÈÉ ÏÅÐÈÎÄ \"!\" , ÂÅÐÎßÒÍÎ, ÍÅ ÑÎÃËÀÑÓÅÒÑß Ñ ÊÎÐÐÅÊÒÈÐÓÅÌÛÌ Â ÄÀÍÍÛÉ ÌÎÌÅÍÒ ÐßÄÎÌ. ÐÅÊÎÌÅÍÄÓÅÒÑß UPDATEFROM ÍÀ×ÈÍÀß Ñ ÄÀÒÛ: \"!\".");
							break;

						case 15:
							strcpy(ret, "Â ÖÅËÅÂÎÌ ÐßÄÅ ÏÐÈÑÓÒÑÒÂÓÞÒ ÎÒÐÈÖÀÒÅËÜÍÛÅ ÇÍÀ×ÅÍÈß");
							break;

						case 16:
							strcpy(ret, "ÐßÄ ÒÅÊÓÙÈÕ ÏÎÊÀÇÀÒÅËÅÉ ÑÎÄÅÐÆÈÒ ÍÓËÅÂÛÅ È/ÈËÈ ÎÒÐÈÖÀÒÅËÜÍÛÅ ÇÍÀ×ÅÍÈß");
							break;

						case 17:
							strcpy(ret, "ÎÒÑÓÒÑÒÂÓÞÙÈÅ ÄÀÍÍÛÅ Â ÊÎÍÖÅ ÊÎÍÒÐÎËÜÍÎÃÎ ÐßÄÀ. ÐßÄ ÒÅÊÓÙÈÕ ÏÎÊÀÇÀÒÅËÅÉ ÝÊÑÒÐÀÏÎËÈÐÎÂÀÍ.");
							break;

						case 18:
							strcpy(ret, "ERROR 10:   ÍÅÂÎÇÌÎÆÍÎ ÏÐÎ×ÈÒÀÒÜ ÖÅËÅÂÎÉ ÏÎÊÀÇÀÒÅËÜ");
							break;

						case 19:
							strcpy(ret, "ERROR 11:   ÍÅÂÎÇÌÎÆÍÎ ÏÐÎ×ÈÒÀÒÜ ÒÅÊÓÙÈÉ ÏÎÊÀÇÀÒÅËÜ. ÏÎÆÀËÓÉÑÒÀ, ÏÐÎÂÅÐÜÒÅ ÄÀÒÛ È ÏÅÐÈÎÄÈ×ÍÎÑÒÜ");
							break;

						case 20:
							strcpy(ret, "ERROR 10:   ÍÅÂÎÇÌÎÆÍÎ ÏÐÎ×ÈÒÀÒÜ ÊÎÍÒÐÎËÜÍÛÉ ÏÎÊÀÇÀÒÅËÜ");
							break;

						case 21:
							strcpy(ret, "ÈÑÏÎËÜÇÎÂÀÍ ÒÎËÜÊÎ ÎÄÈÍ ÊÎÍÒÐÎËÜÍÛÉ ÏÎÊÀÇÀÒÅËÜ");
							break;

						default:
							sprintf(ret, "No message for number s%dm%d", setnum, messnum);
							break;

					}

					break;

				default:
					switch (messnum)
					{
						case  1:
							strcpy(ret, "ERROR 1: INVALID VALUE \"!\" FOR KEYWORD: \"!\"");
							break;

						case  2:
							strcpy(ret, "ERROR 2: CONFLICT BETWEEN BENCHFREQ AND FREQ.");
							break;

						case  3:
							strcpy(ret, "ERROR 3: SPECIFIED TIME SPAN NOT LONG ENOUGH.");
							break;

						case  4:
							strcpy(ret, "ERROR 4: CONFLICT BETWEEN UPDATEFROM AND FROM OR TO DATES.");
							break;

						case  5:
							strcpy(ret, "ERROR 5: DISTRIBUTORID AND TARGETID ARE NOT ALLOWED TO\nHAVE THE SAME NAME.  PLEASE RE-ENTER TARGETID.");
							break;

						case  6:
							strcpy(ret, "ERROR 6: UNKNOWN KEYWORD=VALUE PAIR : \"!\" = \"!\".");
							break;

						case  7:
							strcpy(ret, "ERROR 7: SPECIFIED FISCALLAG MUST BE SMALLER THAN FREQ.");
							break;

						case  8:
							strcpy(ret, "ERROR 8: CANNOT HAVE MEAN=Y AND STOCK=Y AT THE SAME TIME.");
							break;

						case  9:
							strcpy(ret, "ERROR 9: MISSING VALUE FOR KEYWORD: \"!\".");
							break;

						case 10:
							strcpy(ret, "CANNOT MAKE SAME DISTRIBUTOR EQUAL TO TWO DIFFERENT BENCHMARKS. FIRST BENCHMARK IS LINK POINT VALUE. OTHER VALUE DROPPED.");
							break;

						case 11:
							strcpy(ret, "IT IS LIKELY THAT A DISCONTINUITY IN MOVEMENT IS PRESENT BETWEEN THE ADJUSTED SERIES CALCULATED HERE AND THAT WHICH PRECEDES ON THE BASE. RECOMMENDED FROM DATE: \"!\".");
							break;

						case 12:
							strcpy(ret, "IT IS LIKELY THAT THE BENCHMARK COVERING PERIOD \"!\" IS NOT SATISFIED BY THE ADJUSTED SERIES CALCULATED HERE. RECOMMENDED FROM DATE: \"!\".");
							break;

						case 13:
							strcpy(ret, "IT IS LIKELY THAT A DISCONTINUITY IN MOVEMENT IS PRESENT BETWEEN THE UPDATED ADJUSTED SERIES CALCULATED HERE AND THAT WHICH PRECEDES ON THE BASE. RECOMMENDED UPDATEFROM DATE: \"!\".");
							break;

						case 14:
							strcpy(ret, "IT IS LIKELY THAT THE BENCHMARK COVERING PERIOD \"!\" IS NOT SATISFIED BY THE UPDATED ADJUSTED SERIES CALCULATED HERE. RECOMMENDED UPDATEFROM DATE: \"!\".");
							break;

						case 15:
							strcpy(ret, "TARGET SERIES HAS SOME NEGATIVE VALUES");
							break;

						case 16:
							strcpy(ret, "DISTRIBUTOR SERIES HAS SOME VALUES SMALLER OR EQUAL TO 0");
							break;

						case 17:
							strcpy(ret, "MISSING DATA AT THE END OF BENCHMARK SERIES. DISTRIBUTOR EXTRAPOLATED.");
							break;

						case 18:
							strcpy(ret, "ERROR 10:   COULD NOT READ TARGET");
							break;

						case 19:
							strcpy(ret, "ERROR 11:   COULD NOT READ DISTRIBUTOR PLEASE CHECK THE DATES AND THE FREQUENCY");
							break;

						case 20:
							strcpy(ret, "ERROR 10:   COULD NOT READ BENCHMARK");
							break;

						case 21:
							strcpy(ret, "ONLY ONE BENCHMARK USED");
							break;

						default:
							sprintf(ret, "No message for number s%dm%d", setnum, messnum);
							break;

					}

					break;
			}

			break;

		case 3:
			switch (langnum)
			{
				default:
					switch (messnum)
					{
						case  1:
							strcpy(ret, "BENCHFREQ=  1,4                                      1\nFREQ     =  4,12                                     4\nFROM     =  YY, YYQQ or YYMM                         NONE\nTO       =  YY, YYQQ or YYMM                         NONE\nFISCALLAG=  -FREQ < FISCALLAG < FREQ                 0\n");
							break;

						case  2:
							strcpy(ret, "LINKED    = Y-YES, N-NO                              N\nROUND     = Y-YES, N-NO                              N\nDECS      = 0, 1, 2, 3, 4, 5                         1\nPROP      = Y-YES, N-NO      PROPORTIONAL ADJUST.    Y\nFIRST     = Y-YES, N-NO      FIRST DIFF   ADJUST.    Y\nUPDATE    = Y-YES, N-NO                              N\nUPDATEFROM= YY, YYQQ, YYMM                           FROM date\nMEAN      = Y-YES, N-NO      INDEX SERIES            N\nSTOCK     = Y-YES, N-NO      STOCK SERIES            N\n");
							break;

						case  3:
							strcpy(ret, "ARATES =    Y-YES, N-NO                              N\nDIFF   =    Y-YES, N-NO      FIRST DIFFERENCES       N\nFACT   =    Y-YES, N-NO      ADJUSTMENT FACTORS      N\nGR     =    Y-YES, N-NO      GROWTH RATE             N\nLAG    =    1,FREQ                                   1 OR FREQ\n");
							break;

						case  4:
							strcpy(ret, "BENCHID       = VALID SERIES NAME                    NONE\nDISTRIBUTORID = VALID SERIES NAME                    NONE\nTARGETID      = VALID SERIES NAME                    NONE\n");
							break;

						default:
							sprintf(ret, "No message for number s%dm%d", setnum, messnum);
							break;
					}

					break;
			}

			break;

		case 4:
			switch (langnum)
			{
				case 1:
					switch (messnum)
					{
						case  1:
							strcpy(ret, "ENTREZ LES INFORMATIONS SUR LES SERIES\n\nMOT CLEF    VALEURS                                  VALEURS DE DEFAULT");
							break;

						case  2:
							strcpy(ret, "ENTREZ LES OPTIONS DE L'ALGORITHME\n\nMOT CLEF    VALEURS                                  VALEURS DE DEFAULT");
							break;

						case  3:
							strcpy(ret, "CHOISISSEZ LES RAPPORTS\n\nMOT CLEF    VALEURS                                  VALEURS DE DEFAULT");
							break;

						case  4:
							strcpy(ret, "ENTREZ LE NOM DES SERIES A TRAITER");
							break;

						case  5:
							strcpy(ret, "ENTREZ \"MOT-CLE=VALEUR,....\" OU \"E\" POUR TERMINER");
							break;

						case  6:
							strcpy(ret, "SERIE NON AJUSTEE: \"!\"\n\n                TOTAL               I            II           III            IV");
							break;

						case  7:
							strcpy(ret, "AJUSTEMENT DE JALONS UTILISANT LA TECHNIQUE DE MINIMISATION QUADRATIQUE: VERSION DU 28 NOV 90\n\n\n");
							break;

						case  8:
							strcpy(ret, "RAPPORT RELATIF AUX PERIODES FISCALES:\n\nPERIODE           VALEURS CORRESPONDANTES DES                  ECARTS\n                SERIES  AJUSTEE     NON-AJUSTEE  PROPORTIONNELS   ET     ADDITIF");
							break;

						case  9:
							strcpy(ret, "SERIE AJUSTEE: \"!\"\n\n                TOTAL               I            II           III            IV");
							break;

						case 10:
							strcpy(ret, "FACTEUR D'AJUSTEMENTS:\n\n                                    I            II           III            IV");
							break;

						case 11:
							strcpy(ret, "PREMIERES DIFFERENCES DE LA");
							break;

						case 12:
							strcpy(ret, "FACTEUR DE CROISSANCE DE LA");
							break;

						case 13:
							strcpy(ret, "i INDIQUE UNE PERIODE INCOMPLETE.");
							break;

						case 14:
							strcpy(ret, "i INDIQUE UNE PERIODE INCOMPLETE.\nb INDIQUE QU'IL N'Y A PAS DE JALON COUVRANT LA PERIODE.");
							break;

						case 15:
							strcpy(ret, "**** OPTIONS COURANTES ****");
							break;

						case 16:
							strcpy(ret, "RAPPORT RELATIF AUX PERIODES CIVILES:\n\nPERIODE           VALEURS CORRESPONDANTES DES                  ECARTS\n                SERIES  AJUSTEE     NON-AJUSTEE  PROPORTIONNELS   ET     ADDITIF");
							break;

						case 17:
							strcpy(ret, "SEQUENCE NUMERO: \"!\"");
							break;

						case 18:
							strcpy(ret, "SERIE NON AJUSTEE: \"!\"\n\n                                    I            II           III            IV");
							break;

						case 19:
							strcpy(ret, "SERIE AJUSTEE: \"!\"\n\n                                    I            II           III            IV");
							break;

						default:
							sprintf(ret, "No message for number s%dm%d", setnum, messnum);
							break;
					}

					break;

				case 2:
					switch (messnum)
					{
						case  1:
							strcpy(ret, "ÓÏÐÀÂËßÞÙÀß ÈÍÔÎÐÌÀÖÈß\n\nÊËÞ×ÅÂÛÅ ÑËÎÂÀ    ÇÍÀ×ÅÍÈß                                   ÏÎ ÓÌÎË×ÀÍÈÞ");
							break;

						case  2:
							strcpy(ret, "ÀËÃÎÐÈÒÌ ÊÂÀÄÐÀÒÍÎÉ ÌÈÍÈÌÈÇÀÖÈÈ\n\nÊËÞ×ÅÂÛÅ ÑËÎÂÀ    ÇÍÀ×ÅÍÈß                                   ÏÎ ÓÌÎË×ÀÍÈÞ");
							break;

						case  3:
							strcpy(ret, "ÎÒ×ÅÒÛ\n\nÊËÞ×ÅÂÛÅ ÑËÎÂÀ    ÇÍÀ×ÅÍÈß                                   ÏÎ ÓÌÎË×ÀÍÈÞ");
							break;

						case  4:
							strcpy(ret, "ÈÄÅÍÒÈÔÈÊÀÒÎÐÛ ÐßÄÎÂ");
							break;

						case  5:
							strcpy(ret, "ÂÂÅÄÈÒÅ \"KEYWORD=VALUE,....\", \"CONTROL\" ÈËÈ \"E\" ÄËß ÂÛÕÎÄÀ ÈÇ ÏÐÎÃÐÀÌÌÛ");
							break;

						case  6:
							strcpy(ret, "ÍÅÑÊÎÐÐÅÊÒÈÐÎÂÀÍÍÛÉ ÐßÄ: \"!\"\n\n                ÂÑÅÃÎ               I            II           III            IV");
							break;

						case  7:
							strcpy(ret, "ÓÂßÇÛÂÀÍÈÅ ÌÅÒÎÄÎÌ ÊÂÀÄÐÀÒÍÎÉ ÌÈÍÈÌÈÇÀÖÈÈ ÂÅÐÑÈß: 28 ÍÎßÁÐß 1990\n\n\n");
							break;

						case  8:
							strcpy(ret, "ÎÒ×ÅÒ ÇÀ ÔÈÍÀÍÑÎÂÛÉ ÏÅÐÈÎÄ:\n\nÏÅÐÈÎÄ                  ÇÍÀ×ÅÍÈß ÐßÄÎÂ               ÐÀÑÕÎÆÄÅÍÈß\n                       ÑÊÎÐÐÅÊÒÈÐ. ÍÅÑÊÎÐÐ.   ÏÐÎÏÎÐÖÈÎÍÀËÜÍÛÅ   È   ÀÄÄÈÒÈÂÍÛÅ");
							break;

						case  9:
							strcpy(ret, "ÑÊÎÐÐÅÊÒÈÐÎÂÀÍÍÛÉ ÐßÄ: \"!\"\n\n                ÂÑÅÃÎ               I            II           III            IV");
							break;

						case 10:
							strcpy(ret, "ÊÎÝÔÔÈÖÈÅÍÒÛ ÏÎÏÐÀÂÎÊ:\n\n                                    I            II           III            IV");
							break;

						case 11:
							strcpy(ret, "ÐÀÑÕÎÆÄÅÍÈÅ ÏÅÐÂÎÃÎ ÏÎÐßÄÊÀ  ");
							break;

						case 12:
							strcpy(ret, "ÏÐÎÖÅÍÒÍÛÅ ÈÇÌÅÍÅÍÈß");
							break;

						case 13:
							strcpy(ret, "i ÑÂÈÄÅÒÅËÜÑÒÂÓÅÒ Î ÍÅÏÎËÍÎÌ ÈÒÎÃÎÂÎÌ ÏÎÊÀÇÀÒÅËÅ ÇÀ ÏÅÐÈÎÄ.");
							break;

						case 14:
							strcpy(ret, "i ÑÂÈÄÅÒÅËÜÑÒÂÓÅÒ Î ÍÅÏÎËÍÎÌ ÈÒÎÃÎÂÎÌ ÏÎÊÀÇÀÒÅËÅ ÇÀ ÏÅÐÈÎÄ.\nb ÓÊÀÇÛÂÀÅÒ, ×ÒÎ ÄËß ÄÀÍÍÎÃÎ ÏÅÐÈÎÄÀ ÎÒÑÓÒÑÒÂÓÅÒ ÊÎÍÒÐÎËÜÍÛÉ ÏÎÊÀÇÀÒÅËÜ.");
							break;

						case 15:
							strcpy(ret, "**** ÒÅÊÓÙÈÅ ÓÑÒÀÍÎÂÊÈ ****");
							break;

						case 16:
							strcpy(ret, "ÎÒ×ÅÒ ÇÀ ÊÀËÅÍÄÀÐÍÛÉ ÏÅÐÈÎÄ:\n\nÏÅÐÈÎÄ                  ÇÍÀ×ÅÍÈß ÐßÄÎÂ               ÐÀÑÕÎÆÄÅÍÈß\n                       ÑÊÎÐÐÅÊÒÈÐ. ÍÅÑÊÎÐÐ.   ÏÐÎÏÎÐÖÈÎÍÀËÜÍÛÅ   È   ÀÄÄÈÒÈÂÍÛÅ");
							break;

						case 17:
							strcpy(ret, "ÍÎÌÅÐ ÈÒÅÐÀÖÈÈ: \"!\"");
							break;

						case 18:
							strcpy(ret, "ÍÅÑÊÎÐÐÅÊÒÈÐÎÂÀÍÍÛÉ ÐßÄ: \"!\"\n\n                                    I            II           III            IV");
							break;

						case 19:
							strcpy(ret, "ÑÊÎÐÐÅÊÒÈÐÎÂÀÍÍÛÉ ÐßÄ: \"!\"\n\n                                    I            II           III            IV");
							break;

						default:
							sprintf(ret, "No message for number s%dm%d", setnum, messnum);
							break;
					}

					break;

				default:
					switch (messnum)
					{
						case  1:
							strcpy(ret, "CONTROL INFORMATION\n\nKEYWORDS    VALUES                                   DEFAULT");
							break; 

						case  2:
							strcpy(ret, "QUADRATIC MINIMIZATION ALGORITHM\n\nKEYWORDS    VALUES                                   DEFAULT");
							break; 

						case  3:
							strcpy(ret, "REPORTS\n\nKEYWORDS    VALUES                                   DEFAULT");
							break; 

						case  4:
							strcpy(ret, "SERIES IDENTIFICATORS");
							break; 

						case  5:
							strcpy(ret, "ENTER \"KEYWORD=VALUE,....\", \"CONTROL\" OR \"E\" TO EXIT PROGRAM");
							break; 

						case  6:
							strcpy(ret, "UNADJUSTED SERIES: \"!\"\n\n                TOTAL               I            II           III            IV");
							break; 

						case  7:
							strcpy(ret, "BENCHMARK ADJUSTMENT USING QUADRATIC MINIMIZATION TECHNIQUE VERSION: 28 NOV. 90\n\n\n");
							break; 

						case  8:
							strcpy(ret, "FISCAL PERIOD REPORT:\n\nPERIOD                  VALUES FOR THE SERIES               DISCREPANCIES\n                       ADJUSTED      UNADJUSTED   PROPORTIONAL   AND   ADDITIVE");
							break; 

						case  9:
							strcpy(ret, "ADJUSTED SERIES: \"!\"\n\n                TOTAL               I            II           III            IV");
							break; 

						case 10:
							strcpy(ret, "ADJUSTMENT FACTORS:\n\n                                    I            II           III            IV");
							break; 

						case 11:
							strcpy(ret, "FIRST DIFFERENCES OF");
							break; 

						case 12:
							strcpy(ret, "PERCENTAGE CHANGES OF");
							break; 

						case 13:
							strcpy(ret, "i INDICATES AN INCOMPLETE PERIOD TOTAL.");
							break; 

						case 14:
							strcpy(ret, "i INDICATES AN INCOMPLETE PERIOD TOTAL.\nb INDICATES THAT NO BENCHMARK COVERS THIS PERIOD.");
							break; 

						case 15:
							strcpy(ret, "**** CURRENT SETTINGS ****");
							break; 

						case 16:
							strcpy(ret, "CALENDAR PERIOD REPORT:\n\nPERIOD                  VALUES FOR THE SERIES               DISCREPANCIES\n                       ADJUSTED      UNADJUSTED   PROPORTIONAL   AND   ADDITIVE");
							break; 

						case 17:
							strcpy(ret, "RUN NUMBER: \"!\"");
							break; 

						case 18:
							strcpy(ret, "UNADJUSTED SERIES: \"!\"\n\n                                    I            II           III            IV");
							break; 

						case 19:
							strcpy(ret, "ADJUSTED SERIES: \"!\"\n\n                                    I            II           III            IV");
							break; 

						default:
							sprintf(ret, "No message for number s%dm%d", setnum, messnum);
							break;
					}

					break;
			}

			break;

		default:
			sprintf(ret, "No messages for set s%d", setnum);
			break;

	}

	return(ret);
}



void print_default(double *dist, double *trget, char from[], int freq,
	int benchfreq, int nbpoints, int ndecs, int div, char stock, char *prnt);
void print_fisc(double *dist, double *trget, int *tau, int *kappa, int nbpoints,
	int nbbench, int ndecs, int freq, int benchfreq, char from[],
	int div, char stock, char *prnt);
void cal_print(double distsum, double trgetsum, char *format,
	char date[], int ndecs, int freq, int dateinc);
double sumf(double *in, int nb);
void printem(int ndecs, double trgets, double dists, double fac, double dif);
void prnt_data(char start[], int nbpoints, int freq, int nbdecs,
	double *series, char arates, char printsum);

extern void add_date(char *date, int freq, int val);
extern FILE *tables;




/**********
 *
 * void print_default(double *dist, double *trget, char from[], int freq,
 *                    int benchfreq, int nbpoints, int ndecs, int div)
 *
 * Routine to print the default report.
 * It's done in three parts
 * - Print the total of the first incomplete year or quarter.
 * - Print the totals of complete years or quarters
 * - Print the remaining year or quarter
 *
 * The div parameter is used to represent by which number we should divide
 * the sum of target.  We need to do this when the mean=y option has
 * been choosen, because all the resulting numbers have been multiplied
 * by the frequency.  We need to divide by the frequency to produce
 * result.  When mean=n then the div parameter will be one.
 *
 * Only numbers are printed in this routine.  Titles and legend are
 * printed by the calling function (print_report).
 *
 **********/

void print_default(double *dist, double *trget, char from[], int freq,
	int benchfreq, int nbpoints, int ndecs, int div, char stock, char *prnt)
{
	double distsum;
	double trgetsum;
	char date[7];
	int year;
	int per;
	int nbcal;
	int dateinc;
	int inc;
	int lap;
	int stockinc;

	*prnt = NO;
	dateinc = 0;
	stockinc = 0;
	strcpy(date, from);
	lap = 0;
	nbcal = 0;

	if (benchfreq == 1)
		inc = (freq == 12 ? 12 : 4);
	else
		inc = 3;

	per = atoi(from+4);
	year = (atoi(from) - per) / 100;
	lap = (freq - per + 1) % inc;

	if (stock)
	{
		if (lap)
			add_date(date, freq, lap-1);
		else
			add_date(date, freq, inc-1);
	}

	/**********
	* first incomplete period.
	**********/
	if (lap)
	{
		if (stock)
			stockinc = lap - 1;
		else
			dateinc = lap - 1;

		distsum = sumf(&dist[nbcal+stockinc], lap-stockinc) / div;
		trgetsum = sumf(&trget[nbcal+stockinc], lap-stockinc) / div;
		cal_print(distsum, trgetsum, " %s - %s  i", date, ndecs, freq, dateinc);
		nbcal = lap;

		if (stock)
			add_date(date, freq, inc);
		else
			add_date(date, freq, lap);

		*prnt = YES;
	}

	/**********
	* complete periods.
	**********/
	if (stock)
		stockinc = inc - 1;
	else
		dateinc = inc - 1;

	while (nbcal + inc <= nbpoints)
	{
		distsum = sumf(&dist[nbcal+stockinc], inc-stockinc) / div;
		trgetsum = sumf(&trget[nbcal+stockinc], inc-stockinc) / div;
		cal_print(distsum, trgetsum, " %s - %s   ", date, ndecs, freq, dateinc);
		nbcal += inc;
		add_date(date, freq, inc);
	}

	/**********
	* last incomplete period.
	**********/
	if (nbcal < nbpoints  && !stock)
	{
		inc = nbpoints - nbcal;
		dateinc = inc - 1;
		distsum = sumf(&dist[nbcal], inc) / div;
		trgetsum = sumf(&trget[nbcal], inc) / div;
		cal_print(distsum, trgetsum, " %s - %s  i", date, ndecs, freq, dateinc);
		nbcal += inc;
		*prnt = YES;
	}
}



/**********
 *
 * void print_fisc(double *dist, double *trget, int *tau, int *kappa, int nbpoints
 *	int nbbench, int ndecs, int freq, int benchfreq, char from[],
 *	int div, char stock)
 *
 * Routine to print the fiscal report.
 * It's done in five parts
 * - Print the total of the first incomplete year or quarter.
 * - Print totals of periods for which there is no benchmarks
 * - Print the totals of complete years or quarters
 * - Print totals of periods for witch there is no benchmarks
 * - Print the remaining year or quarter
 *
 * Only numbers are printed in this routine.  Titles and legend are
 * printed by the calling function (print_report).
 * Tau and Kappa are the reference period of the benchmark over the
 * distributor.  Tau is start and kappa is end.
 *
 **********/

void print_fisc(double *dist, double *trget, int *tau, int *kappa, int nbpoints,
	int nbbench, int ndecs, int freq, int benchfreq, char from[],
	int div, char stock, char *prnt)
{
	double distsum;
	double trgetsum;
	char date[7];
	int start;
	int i, stockinc;
	int inc, dateinc;
	int nbcal;
	int lap;
	int nbref;

	*prnt = NO;
	dateinc = 0;
	strcpy(date, from);
	lap = 0;
	nbcal = 0;

	if (benchfreq == 1)
		inc = (freq == 12 ? 12 : 4);
	else
		inc = 3;

	lap = (tau[0] - 1) % inc;

	if (lap && !stock)
	{
		dateinc = lap - 1;
		distsum = sumf(&dist[0], lap) / div;
		trgetsum = sumf(&trget[0], lap) / div;
		cal_print(distsum, trgetsum, " %s - %s bi", date, ndecs, freq, dateinc);
		*prnt = YES;
	}

	add_date(date, freq, lap);

	if (!stock)
		dateinc = inc - 1;

	while (lap < tau[0] - 1)
	{
		stockinc = 0;
		if (stock)
			stockinc = inc - 1;

		distsum = sumf(&dist[lap], inc-stockinc) / div;
		trgetsum = sumf(&trget[lap], inc-stockinc) / div;
		cal_print(distsum, trgetsum, " %s - %s b ", date, ndecs, freq, dateinc);
		lap += inc;
		add_date(date, freq, inc);
		*prnt = YES;
	}

	for (i = 0; i < nbbench; i++)
	{
		start = tau[i] - 1;
		nbref = kappa[i] - start;
		distsum = sumf(&dist[start], nbref) / div;
		trgetsum = sumf(&trget[start], nbref) / div;
		cal_print(distsum, trgetsum, " %s - %s   ", date, ndecs, freq, dateinc);
		add_date(date, freq, inc);
	}

	lap = kappa[nbbench - 1];

	while (lap + inc <= nbpoints)
	{
		distsum = sumf(&dist[lap], inc) / div;
		trgetsum = sumf(&trget[lap], inc) / div;
		cal_print(distsum, trgetsum, " %s - %s b ", date, ndecs, freq, dateinc);
		lap += inc;
		add_date(date, freq, inc);
		*prnt = YES;
	}


	if (lap < nbpoints && !stock)
	{
		inc = nbpoints - lap;
		dateinc = inc - 1;
		distsum = sumf(&dist[lap], inc) / div;
		trgetsum = sumf(&trget[lap], inc) / div;
		cal_print(distsum, trgetsum, " %s - %s bi", date, ndecs, freq, dateinc);
		*prnt = YES;
	}

	if (stock)
	{
		lap = kappa[nbbench-1] + inc;
		if (lap < nbpoints)
		{
			distsum = dist[lap-1] / div;
			trgetsum = trget[lap-1] / div;
			cal_print(distsum, trgetsum, " %s - %s b ", date, ndecs, freq, dateinc);
			*prnt = YES;
		}
	}
}



/**********
 *
 * void cal_print(double distsum, double trgetsum, char *format,
 *	char date[], int ndecs, int freq, int dateinc)
 *
 * calculates adjustment factor, adjusted differences and prints
 * the numbers.
 *
 **********/

void cal_print(double distsum, double trgetsum, char *format,
	char date[], int ndecs, int freq, int dateinc)
{
	double adj_fac, adj_dif;
	char date2[7];

	strcpy(date2, date);
	add_date(date2, freq, dateinc);

	if (trgetsum && distsum)
		adj_fac = trgetsum / distsum;
	else
		adj_fac = 0;

	adj_dif = trgetsum - distsum;
	fprintf(tables, format, date, date2);
	printem(ndecs, trgetsum, distsum, adj_fac, adj_dif);
}



/**********
 *
 * double sumf(double *in, int nb)
 *
 * sums up the value of an array of doubles.
 * the result is put into a double for greater accuracy.
 *
 **********/

double sumf(double *in, int nb)
{
	double total;
	int   i;

	total = 0.0;
	for (i = 0; i < nb; i++)
		total += in[i];

	return(total);
}



/**********
 *
 * void printem(int ndecs, double trgets, double dists, double fac, double dif)
 *
 * print four doubles on a line.
 *
 **********/

void printem(int ndecs, double trgets, double dists, double fac, double dif)
{
	fprintf(tables, "%16.*f", ndecs, trgets);
	fprintf(tables, "%16.*f", ndecs, dists);
	fprintf(tables, "%15.4f", fac);
	fprintf(tables, "%17.4f\n", dif);
}



/**********
 *
 * void prnt_data(char start[], int nbpoints, int freq, int nbdecs,
 *                double *series, char arates, char printsum);
 *
 * This procedure prints a series on a line by line format.
 * It can print a series any frequency (1,4,12).
 * It prints the series with the numbers of decimals passed in
 * parameter.  The parameter start gives the start date of which the
 * series represent. The parameters nbpoints gives the number of
 * datapoints to be printed.
 *
 **********/

void prnt_data(char start[], int nbpoints, int freq, int nbdecs,
	double *series, char arates, char printsum)
{
	int per;
	int year;
	int num;
	int num1;
	int i;
	int j;
	int mult;

	mult = (arates ? freq : 1);
	per = atoi(start+4);
	year = (atoi(start) - per) / 100;
	num = 0;

	switch (freq)
	{
		case 1:
			for (i = 0; i < nbpoints; i++)
				fprintf(tables, " %4.4d%16.*f\n", year++, nbdecs, series[num++]);

			break;

		case 4:
			/**********
			* first incomplete period
			**********/

			if (per != 1)
			{
				fprintf(tables, "\n %4.4d%16s  ", year++, " ");

				for (i = 0; i < 4; i++)
				{
					if (i < per -1 || num >= nbpoints)
						fprintf(tables, "%14s", "-"); 
					else
						fprintf(tables, "%14.*f", nbdecs, series[num++]*(double)mult);
				}

				fprintf(tables, "\n");
			}

			/**********
			* complete periods
			**********/

			while (num + 4 <= nbpoints)
			{
				fprintf(tables, "\n %4.4d", year++);

				if (printsum)
					fprintf(tables,"%16.*f  ", nbdecs, sumf(&series[num],4));
				else
					fprintf(tables,"%16s  ", " ");

				for (i = 0; i < 4; i ++)
					fprintf(tables, "%14.*f", nbdecs, series[num++]*(double)mult);

				fprintf(tables, "\n");
			}

			/**********
			* last incomplete period
			**********/

			if (num < nbpoints)
			{
				fprintf(tables, "\n %4.4d%16s  ", year++, " ");
				for (i = 0; i < 4; i++)
				{
					if (num < nbpoints)
						fprintf(tables, "%14.*f", nbdecs, series[num++]*(double)mult);
					else
						fprintf(tables, "%14s", "-");
				}

				fprintf(tables, "\n");
			}

			break;

		case 12:
			/**********
			* first incomplete period
			**********/

			if (per != 1)
			{
				fprintf(tables, "\n %4.4d%16s  ", year++, " ");

				for (i = 0; i < 3; i++)
				{
					num1 = i;

					if (i != 0)
						fprintf(tables, " %4s%16s  ", " ", " ");

					for (j = 0; j < 4; j++)
					{
						if (num1 < per - 1 || num1 - per + 2 > nbpoints)
							fprintf(tables, "%14s", "-");
						else
						{
							fprintf(tables, "%14.*f", nbdecs,
							series[num1 - per + 1]*(double)mult);
							num++;
						}

						num1 += 3;
					}

					fprintf(tables, "\n");
				}
			}

			/**********
			* complete periods
			**********/

			while (num + 12 <= nbpoints)
			{
				fprintf(tables, "\n %4.4d", year++);

				if (printsum)
					fprintf(tables,"%16.*f  ", nbdecs, sumf(&series[num],12));
				else
					fprintf(tables,"%16s  ", " ");

				for (i = 0; i < 3; i++)
				{
					if (i != 0)
						fprintf(tables, " %4s%16s  ", " ", " ");

					num1 = i;

					for (j = 0; j < 4; j++)
					{
						fprintf(tables, "%14.*f", nbdecs, series[num+num1]*(double)mult);
						num1 += 3;
					}

					fprintf(tables, "\n");
				}

				num += 12;
			}

			/**********
			* last incomplete period
			**********/

			if (num < nbpoints)
			{
				fprintf(tables, "\n %4.4d%16s  ", year++, " ");

				for (i = 0; i < 3; i++)
				{
					num1 = i;

					if (i != 0)
						fprintf(tables, " %4s%16s  ", " ", " ");

					for (j = 0; j < 4; j++)
					{
						if (num + num1 < nbpoints)
							fprintf(tables, "%14.*f", nbdecs, series[num + num1]*(double)mult);
						else
							fprintf(tables, "%14s", "-");

						num1 += 3;
					}

					fprintf(tables, "\n");
				}
			}

			break;
	}
}



int distribround(double *in, double sum, int nvalue, int ndec, double *out);
void roundd(int ndec, int npts, double *in, double *out);
extern void shellsort(int *n, char *type, double *x, int *nlmnt);



/**********
 *
 * int distribround(double *in, double sum, int nvalue, int ndec, double *out)
 *
 *     PURPOSE:
 *
 *        TO ROUND A VECTOR OF VALUES AND DISTRIBUTE THE VALUES
 *        SO THAT THE SUM OF THE ROUNDED VALUES MATCHES A SUPPLIED
 *        VALUE AFTER ROUNDING.
 *
 *     INPUTS:
 *
 *        IN      - VECTOR OF VALUES TO BE ROUNDED & DISTRIBUTED.
 *        SUM     - VALUE TO WHICH THE SUM OF 'IN' IS TO BE MATCHED.
 *                  THIS VALUE IS NOT NECESSARILY ROUNDED.
 *        NVALUE  - # OF VALUES TO BE ROUNDED & DISTRIBUTED.
 *        NDEC    - # OF DECIMAL PLACES REQUESTED
 *        STCROUNDING   .FALSE. - REGULAR ROUNDING METHOD.
 *                      .TRUE.  - STATISTICS CANADA APPROVED METHOD.
 *
 *     OUTPUTS:
 *
 *        OUT     - VECTOR OF ROUNDED VALUES,CORRESPONDING TO 'IN'
 *                  SUCH THAT SUMMATION OF OUT EQUALS 'SUM' ROUNDED TO
 *                  NDEC DECIMAL PLACES
 *        ERROR   - ERROR CODE
 *                     0 SUCCESSFUL
 *                     1 DABS(SUMMATION OF IN VALUES - SUM) > NVALUES
 *
 *     METHOD:
 *
 *        THE DEVIATION OF EACH ROUNDED VALUE FROM THE ORIGINAL VALUE
 *        IS RANKED IN DESCENDING ORDER.  IF AN ADJUST OF +1(-1) MUST
 *        BE MADE TO THE LAST DIGIT TO BE KEPT THEN IT IS ADDED TO
 *        THE ROUNDED VALUE WITH THE SMALLEST(LARGEST) DEVIATION.
 *        THIS PROCESS IS REPEATED UNTIL ALL OF THE DIFFERENCE BETWEEN
 *        'SUM' AND THE SUM OF THE ROUNDED 'IN' VALUES HAS BEEN
 *        REMOVED.
 *
 *        GIVEN A VECTOR OF N TERMS X(I) AND  A VALUE A TO WHICH
 *        SUMMATION OF X MUST AGREE.  LET ~X(I) BE THE VALUE OF X(I)
 *        ROUNDED TO D DECIMAL PLACES.
 *
 *        BEGIN
 *           FORM ~X(I) FOR ALL I
 *           D = ~A - SUMMATION OF ~X(I)
 *           FORM S(I) = ~X(I) - X(I)
 *           IF D > .5 * 10 ** -DEC THEN
 *              FOR J = 1 TO |D| * D **DEC
 *                  IF D > 0 THEN
 *                     FIND THE NEXT SMALLEST S(I) ,CALL IT S(K)
 *                  ELSE
 *                     FIND THE NEXT LARGEST S(I) CALL IT S(K)
 *                  ENDIF
 *                  ~X(K) = ~X(K) + DSIGN(D) * 10 ** -DEC
 *              ENDFOR
 *           ENDIF
 *        END
 *
 **********/

int distribround(double *in, double sum, int nvalue, int ndec, double *out)
{
	double lshift, rshift;
	double roundsum[1];
	double calsum;
	double diff[300];
	double sumdiff;
	double adjustment;
	double adjcount;
	double a, b;
	int ret_val;
	int seqno[300];
	int adjustcount;
	int i, j, z;

	if (nvalue > 300)
	{
		printf("nvalue > 300  = %d\n", nvalue);
		return(1);
	}

	ret_val = 0;
	a = ndec;
	b = -ndec;
	lshift = pow(10.0, a);
	rshift = pow(10.0, b);
	roundsum[0] = sum;
	roundd(ndec, 1, roundsum, roundsum);
	roundd(ndec, nvalue, in, out);

	calsum = 0;
	for (i = 0; i < nvalue; i++)
	{
		calsum = calsum + out[i];
		diff[i] = out[i] - in[i];
	}
	sumdiff = roundsum[0] - calsum;

	if (fabs(sumdiff) < (0.5 * rshift))
	{
		return(ret_val);
	}

	adjcount = (fabs(sumdiff) + 0.5 * rshift) * lshift;
	if (adjcount > nvalue + 1.0)
	{
		return(1);
	}

	adjustcount = (int) adjcount;
	if (adjustcount > nvalue)
	{
		return(1);
	}

	shellsort(&nvalue, "D", diff, seqno);
	adjustment = (sumdiff < 0 ? -1.0 : 1.0);
	adjustment *= rshift;

	/* DEBUG
	printf("distribround 5.0: after shellsort: sumdiff=%12.4f adjustment=%12.4f\n",sumdiff,adjustment);
	for (z=0; z < nvalue; z++) printf("out[%3d]=%12.4f  diff[%3d]=%14.6f  seqno[%3d]=%3d\n",z,out[z],z,diff[z],z,seqno[z]);
	*/

	if (sumdiff > 0)
	{
		j = nvalue - 1;
		for (i = 0; i < adjustcount; i++)
		{
			out[seqno[j]] = out[seqno[j]] + adjustment;
			j--;
		}
	}
	else
	{
		for (i = 0; i < adjustcount; i++)
		{
			out[seqno[i]] = out[seqno[i]] + adjustment;
		}
	}

	/* DEBUG
	printf("distribround 6.0:\n");
	for (z=0; z < nvalue; z++) printf("out[%3d]=%12.4f\n",z,out[z]);
	*/

	return(ret_val);
}



/**********
 *
 * void roundd(int ndec, int npts, double *in, double *out)
 *
 *     FUNCTION:
 *        TO ROUND A VECTOR OF DATA VALUES TO A GIVEN
 *        # OF DECIMAL PLACES.
 *
 *     INPUTS:
 *        NDEC - # OF DECIMAL PLACES
 *        NPTS - # OF DATA VALUES IN 'IN','OUT'
 *        IN   - VECTORS CONTAINING THE DATA VALUES TO BE
 *               ROUNDED
 *        STCROUNDING - LOGICAL VALUE TO SPECIFY ROUNDING METHOD.
 *                      .FALSE. - REGULAR ROUNDING METHOD.
 *                      .TRUE.  - STATISTICS CANADA APPROVED METHOD.
 *
 *     OUTPUT:
 *        OUT  - VECTOR CONTAINING THE ROUNDED DATA VALUES
 *
 *
 **********/

void roundd(int ndec, int npts, double *in, double *out)
{
	double realnum;
	double intnum;
	double dfactor;
	double lshift;
	double rshift;
	char cfactor[128];
	int i;


	lshift = pow(10.0, (double)ndec);
	rshift = pow(10.0, (double)-ndec);

	for (i = 0; i < npts; i++)
	{
		realnum = in[i] * lshift;
		intnum = floor(realnum);
		if (realnum < 0.0)
			intnum++;

		dfactor = fabs(realnum - intnum);
		sprintf(cfactor, "%lf", dfactor);
		dfactor = atof(cfactor);

		if (dfactor >= 0.5)
		{
			if (realnum < 0)
				intnum--;
			else
				intnum++;
		}

		out[i] = intnum * rshift;
	}
}



/*
C
C THIS SUBROUTINE PERFORMS A SHELL SORT OF A DATA SERIES
C where nlmnt contains the sorted index of the data in
C the x array.
C IF TYPE EQUALS "A", ROUTINE SORTS IN ASCENDING ORDER
C IF TYPE EQUALS "D", ROUTINE SORTS IN DESCENDING ORDER
C
*/

void shellsort(int *n, char *type, double *x, int *nlmnt)
{
	int igap, imax, iex, iplusg, temp, i;

	for (i=0; i<*n; i++)
	{
		nlmnt[i]=i;
	}

	igap=*n;

	if (strcmp(type, "A") == 0)
	{
		while (igap > 1)
		{
			igap = igap/2;
			imax = *n-igap;
			do
			{
				iex=0;
				for (i=0; i<imax; i++)
				{
					iplusg = i + igap;
					if (x[nlmnt[i]] > x[nlmnt[iplusg]])
					{
						temp = nlmnt[i];
						nlmnt[i] = nlmnt[iplusg];
						nlmnt[iplusg] = temp;
						iex++;
					}
				}
			} while (iex > 0);
		}
	}
	else
	{
		while (igap > 1)
		{
			igap = igap/2;
			imax = *n-igap;
			do
			{
				iex=0;
				for (i=0; i<imax; i++)
				{
					iplusg = i + igap;
					if (x[nlmnt[i]] < x[nlmnt[iplusg]])
					{
						temp = nlmnt[i];
						nlmnt[i] = nlmnt[iplusg];
						nlmnt[iplusg] = temp;
						iex++;
					}
				}
			} while (iex > 0);
		}
	}
}

