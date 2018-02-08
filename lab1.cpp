#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <cmath>

using namespace std;

int k;
float threshold;
int iterations;
int workers;
char* input;

class Point {
public:
	vector<float> coordinates;
	int label;

	Point(vector<float>& coordinates) {
		int size = coordinates.size();
		for(int i = 0; i < size; i++) {
			this->coordinates.push_back(coordinates[i]);
		}
	}

	float getDimensionsCount() {
		return coordinates.size();
	}
};

class Centroid {
public:
	int label;
	vector<Point> points;
	vector<float> coordinates;

	Centroid(int label, Point point) {
		this->label = label;
		int count = point.getDimensionsCount();
		for(int i = 0; i < count; i++) {
			coordinates.push_back(point.coordinates[i]);
		}
		points.push_back(point);
	}
};

vector<Centroid> randomCentroids(vector<Point> points, int k) {
	vector<Centroid> retval;

	for(int i = 0; i < k; i++) {
		while (1) {
			int index = rand() % points.size();
			if (i > 0) {
				if (retval[i-1].coordinates == points[i].coordinates) {
					break;
				}
			}

			Centroid* centroid = new Centroid(i, points[index]);
			retval.push_back(*centroid);
			break;
		}
	}

	return retval;
}


void findNearestCentroids(vector<Point> points, vector<Centroid> centroids) {
	int total_points = points.size();
	int total_coordinates = points[0].getDimensionsCount();

	for (int i = 0; i < total_points; i++) {
		Point curr_point = points[i];

		float sum_diff = 0;
		float min;

		for (int j = 0; j < total_coordinates; j++) {
			sum_diff += pow(centroids[0].coordinates[j] - curr_point.coordinates[j], 2.0);
		}

		min = sqrt(sum_diff);

		for (int j = 1; j < k; j++) {
			float dist;
			float sum = 0;

			for(int m = 0; m < total_coordinates; m++) {
				sum += pow(centroids[j].coordinates[m] - curr_point.coordinates[m], 2.0);
			}

			dist = sqrt(sum);

			if(min > dist) {
				min = dist;
				curr_point.label = i;
			}
		}
	}
}

vector<Centroid> averageLabeledCentroids(vector<Point> points, vector<Centroid> centroids) {
	int total_points = points.size();
	vector<Centroid> retval;

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < total_points; j++) {
			float sum = 0;

			for (int m = 0; m < total_points; m++) {
				sum += centroids[i].points[m].coordinates[j];
			}
			centroids[i].coordinates[j] = sum / total_points;
		}
		retval.push_back(centroids[i]);
	}

	return retval;
}

void kmeans(vector<Point> dataset, int k) {
	vector<Centroid> centroids = randomCentroids(dataset, k);

	int iters = 0;
	vector<Centroid> old_centroids;

	bool done = false;

	while (!done) {
		old_centroids = centroids;
		iters++;

		findNearestCentroids(dataset, centroids);
		centroids = averageLabeledCentroids(dataset, centroids);

		done = iters > iterations;
	}
}

int main (int argc, char **argv) {
	int c;
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
