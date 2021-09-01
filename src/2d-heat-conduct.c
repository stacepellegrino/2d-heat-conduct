#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <getopt.h>
#include <string.h>

void declare_2d_array(double ***array, int rows, int cols);

void free_2d_array(double ***array, int size);

static char* algorithm = NULL;
static int chunk = 0;
static double degrees[4] = {100.0, 100.0, 100.0, 100.0}; // left, right, top, bottom temperatures.
static omp_sched_t sched = 1;
static int nthreads = 1;
static int noprint = 0;
static int sizes[2] = {100, 100}; // sizes x and y
static double tolerance = 0.0001;
static int wtime = 0;
static int verbose = 0;
static const char* version = "2d-heat-conduct 1.0.0";
static const char* usage = "usage: 2d-heat-conduct [-hpvVw] [-a algorithm] [-c chunk-size] [-d left-temp, right-temp, top-temp, bottom-temp]  [-s x y] [-S schedule] [-t tolerance] [-T threads]";

int main(int argc, char **argv)
{
	double tstart, tstop;

	tstart = omp_get_wtime();

	int c;

	int index = 0;
	int dcount = 0;
	int scount = 0;
	char* next = NULL;

	char* schedstr = "static";

	int m; 
	int n;

	int i, j, iter;

	
	void declare_2d_array(double ***array, int rows, int cols) {
		*array = malloc((rows) * sizeof(double *));
		if (*array == NULL)
		{
			fprintf(stderr, "No memory\n");
			exit(1); 
		}

		for (int i = 0; i < rows; i++)
		{
			(*array)[i] = malloc((cols) * sizeof(double));
			if ((*array)[i] == NULL)
			{
				fprintf(stderr, "No memory\n");
				exit(1); 
			}
			for (int j = 0; j < cols; j++) {
				(*array)[i][j] = 30.0;
			}
		}	
	}

	void free_2d_array(double ***array, int size) {
		for (int i = 0; i < size; i++)
		{
			free((*array)[i]);
			(*array)[i] = NULL;
		}
		free(*array);
		*array = NULL;
	}

	while(1)
	{
		static struct option long_options[] =
		{
			{"algorithm", required_argument, 0, 'a'},
			{"chunk",     required_argument, 0, 'c'},
			{"degrees",   required_argument, 0, 'd'},
			{"help",      no_argument,       0, 'h'},
			{"noprint",   no_argument,       0, 'p'},
			{"size",      required_argument, 0, 's'},
			{"schedule",  required_argument, 0, 'S'},
			{"tolerance", required_argument, 0, 't'},
			{"threads",   required_argument, 0, 'T'},
			{"verbose",   no_argument,       0, 'v'},
			{"version",   no_argument,       0, 'V'},
			{"wtime",     no_argument,       0, 'w'},
			{NULL, 0, NULL, 0}
		};

		int option_index = 0;

		c = getopt_long (argc, argv, "a:c:d:hps:S:t:T:vVw", long_options, &option_index);

		// Detect the end of the options.
		if (c == -1)
			break;

		switch(c)
		{
			case 'a':
				if (strcmp(optarg, "jacobi") == 0 || strcmp(optarg, "gauss") == 0) {
					algorithm = optarg;
				}
				else  {
					fprintf(stderr, "error: -a requires [jacobi] or [gauss] as an arg\n");
					exit(1);
				}
				break;

			case 'c':
				if (atoi(optarg) >= 0 && atoi(optarg) <= 256) {
					chunk = atoi(optarg);
				}
				else {
					fprintf(stderr, "error: -c requires a chunk size from 0 to 256 an arg\n");
					exit(1);
				}
				break;

			case 'd':
				index = optind - 1;
				while(index < argc) { 
					next = strdup(argv[index]);
					index++;
					if (atof(next) > 0.000000 && atof(next) <= 1000.0) {
						degrees[dcount++] = atof(next);
                			}
			                else break;

            			}
				if (dcount != 4) {
					fprintf(stderr, "error: -d requires four temperatures (above 0 up to 1000.0 degrees) in decimal as args [left] [right] [top] [bottom]\n");
					exit(1);
				}
				if (argc == optind + 3) {
					optind = index;
				}
				else optind = index - 1;
				break;

			case 'h':
				printf("%s\n", usage);
				exit(0);

			case 'p':
				noprint = 1;
				break;

			case 's':
				index = optind - 1;
				while(index < argc) { 
					next = strdup(argv[index]);
					index++;
					if (atoi(next) > 0 && atoi(next) <= 1000) {
						sizes[scount++] = atoi(next);
                			}
			                else break;

            			}
				if (scount != 2) {
					fprintf(stderr, "error: -s requires two integer sizes (between 1 to 1000000) as args [x y]\n");
					exit(1);
				}
				if (argc == optind + 1) {
					optind = index;
				}
				else optind = index - 1;
				break;

			case 'S':
				if (strcmp(optarg, "static") == 0) {
					sched = 1;
				}
				else if (strcmp(optarg, "dynamic") == 0) {
					sched = 2;
				}
				else if (strcmp(optarg, "guided") == 0) {
					sched = 3;
				}
				else if (strcmp(optarg, "auto") == 0) {
					sched = 4;
				}
				else {
					fprintf(stderr, "error: -S requires [static], [dynamic], [guided] or [auto] schedule as arg\n"); 
					exit(1);
				}
				schedstr = optarg;
				break;

			case 't':
				if (atof(optarg) > 0.00001 && atof(optarg) <= 1.0) {
					tolerance = atof(optarg);
				}
				else {
					fprintf(stderr, "error: -t requires tolerance decimal as arg\n");
					exit(1);
				}
				break;

			case 'T':
				if (atoi(optarg) > 0 && atoi(optarg) <= 256) {
					nthreads = atoi(optarg);
				}
				else {
					fprintf(stderr, "error: -T requires number of threads from 1 to 256 as arg\n");
					exit(1);
				}	
				break;

			case 'v':
				verbose = 1;
				break;

			case 'V':
				printf("%s\n", version);
				exit(0);

			case 'w':
				wtime = 1;
				break;

			case '?':
				fprintf(stderr, "%s\n", usage);
				exit(1);

			default:
				abort();
		}
	} // while end

	if (algorithm == NULL) {
		fprintf(stderr, "error: [-a algorithm] not set\n");
		exit(1);
	}
	

	if (argc > optind || argc == 1) {
		fprintf(stderr, "error: only takes flag options with no further arguments\n\n%s\n", usage);
		exit(1);
	}

	if (verbose) {

		printf("algorithm: %s\n", algorithm);
		printf("chunk: %d\n", chunk);
		printf("degrees: %f, %f, %f, %f\n", degrees[0], degrees[1], degrees[2], degrees[3]);
		printf("noprint: %d\n", noprint);
		printf("sizes: %d, %d\n", sizes[0], sizes[1]);
		printf("schedule: %s\n", schedstr);
		printf("tolerance: %f\n", tolerance);
		printf("threads: %d\n", nthreads);
		printf("verbose: %d\n", verbose);
		printf("wtime: %d\n", wtime);

	}

	omp_set_num_threads(nthreads);
	omp_set_schedule(sched, chunk);

	m = sizes[0];
	n = sizes[1];

	double **t;
	double **tnew;

	declare_2d_array(&t, m + 2, n + 2);
	declare_2d_array(&tnew, m + 2, n + 2);

	double diff, difmax;

	#pragma omp parallel default(shared) private(i, j)
	{
		// boundary conditions
		#pragma omp for nowait
		for (i = 1; i <= m; i++) {
			t[i][0] = degrees[0];
			t[i][n + 1] = degrees[1];
		}

		#pragma omp for nowait
		for (j = 1; j <= n; j++) {
			t[0][j] = degrees[2];
			t[m + 1][j] = degrees[3];
		}
	}

	iter = 0;
	difmax = 1000000.0;

	if (strcmp(algorithm, "jacobi") == 0) {
	
		while (difmax > tolerance) {

			iter++;
			difmax = 0.0;

			#pragma omp parallel default(shared) private(i, j, diff) reduction(max:difmax)
			{
				#pragma omp for nowait collapse(2)
				for (i = 1; i <= m; i++) {
					for (j = 1; j <= n; j++) {
						tnew[i][j] = (t[i - 1][j] + t[i + 1][j] + t[i][j - 1] + t[i][j + 1]) / 4.0;
					}
				}
					
				#pragma omp for nowait collapse(2)
				for (i = 1; i <= m; i++) {
					for (j = 1; j <= n; j++) {
						// update temperature for next iteration
						tnew[i][j] = (t[i - 1][j] + t[i + 1][j] + t[i][j - 1] + t[i][j + 1]) / 4.0;
						diff = fabs(tnew[i][j] - t[i][j]);
						if (diff > difmax) {
							difmax = diff;
						}
						// copy new to old temperatures
						t[i][j] = tnew[i][j];
					}
				}
			}

		}
	}
	else if (strcmp(algorithm, "gauss") == 0) {

		while (difmax > tolerance) {

			iter++;
			difmax = 0.0;

			#pragma omp parallel default(shared) private(i, j, diff) reduction(max:difmax)
			{
				#pragma omp for nowait collapse(2)
				for (i = 1; i <= m; i++) {
					for (j = 1; j <= n; j++) {
						// update temperature for next iteration
						tnew[i][j] = (t[i - 1][j] + t[i + 1][j] + t[i][j - 1] + t[i][j + 1]) / 4.0;
						diff = fabs(tnew[i][j] - t[i][j]);
						if (diff > difmax) {
							difmax = diff;
						}
						// copy new to old temperatures
						t[i][j] = tnew[i][j];
					}
				}
			}

		}
	}


	tstop = omp_get_wtime();

	if (noprint == 0) {

		// print results
		printf("iter = %d  difmax = %9.11lf", iter, difmax);
		for (i = 0; i <= m + 1; i++) {
			printf("\n");
			for (j = 0; j <= n + 1; j++) {
				if (strcmp(algorithm, "jacobi") == 0) {
					printf("%3.5lf ", t[i][j]);
				}
				else if (strcmp(algorithm, "gauss") == 0) {
					printf("%3.5lf ", tnew[i][j]);
				}
			}
		}
		printf("\n");

		printf("Iterations = %d\nMaximum difference = %-5.7lf\n", iter, difmax);
    
	}

	if (wtime) {
		printf("Time taken = %4.7lf seconds\n", (tstop - tstart));
	}

	free_2d_array(&t, m + 2);
	free_2d_array(&tnew, m + 2);

	return 0;
}
