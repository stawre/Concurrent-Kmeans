#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

int main (int argc, char **argv) {
	int c;
	int k;
	float threshold;
	int iterations;
	int workers;
	char* input;
	int option_index = 0;

	static struct option long_options[] = {
		{"input", required_argument, 0, 'i'},
		{"workers", required_argument, 0, 'w'},
		{"threshold", required_argument, 0, 't'},
		{"iterations", required_argument, 0, 'n'},
		{"clusters", required_argument, 0, 'c'},
		{0, 0, 0, 0}
	};

	while ((c = getopt_long(argc, argv, "i:w:t:n:c:", long_options, &option_index)) != -1) {
		switch (c) {
			case 'i':
				input = optarg;
				break;

			case 'w':
				workers = atoi(optarg);
				break;

			case 't':
				threshold = atof(optarg);
				break;

			case 'n':
				iterations = atoi(optarg);
				break;

			case 'c':
				k = atoi(optarg);
				break;

			default:
				abort();
		}
	}

	printf("k: %d, threshold: %f, iterations: %d, workers: %d, input: %s\n", k, threshold, iterations, workers, input);
}