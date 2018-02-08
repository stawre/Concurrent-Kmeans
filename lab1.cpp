#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <vector>

using namespace std;

class Point {
private:
	vector<float> coordinates;
	int label;

public:
	Point(vector<float>& coordinates) {
		int size = coordinates.size();
		for(int i = 0; i < size; i++) {
			this->coordinates.push_back(coordinates[i]);
		}
	}

	float getCoordinate(int index) {
		return coordinates[index];
	}

	float getDimensionsCount() {
		return coordinates.size();
	}
};

class Cluster {
private:
	int label;
	vector<Point> points;
	vector<float> location;

public:
	Cluster(int label, Point point) {
		this->label = label;
		int count = point.getDimensionsCount();
		for(int i = 0; i < count; i++) {
			location.push_back(point.getCoordinate(i));
		}
		points.push_back(point);
	}

	Point getPoint(int index) {
		return points[index];
	}
};

vector<Point> randomCentroids(vector<Point> points, int k) {
	vector<Point> retval;

	for(int i = 0; i < k; i++) {
		while (1) {
			int index = rand() % points.size();
			if (i > 0) {
				if (retval[i-1] == points[i]) {
					break;
				}
			}

			retval.push_back(points[index]);
			Cluster cluster(i, points[index]);
			break;
		}
	}

	return retval;
}


void findNearestCentroids(vector<Point> points, vector<Point> centroids) {
	int total_points = points.size();
	int total_coordinates = points.getDimensionsCount();

	for (int i = 0; i < total; i++) {
		Point curr_point = points[i];

		float sum_diff = 0;
		float min;

		for (int j = 0; j < total_coordinates; j++) {
			sum_diff += pow(clusters[0].location[j] - curr_point.coordinates[j], 2.0);
		}

		min = sqrt(sum);

		for (int j = 1; j < k; j++) {
			float dist;
			float sum = 0;

			for(int m = 0; m < total_coordinates; m++) {
				sum += pow(clusters[j].location[m] - curr_point.coordinates[m], 2.0);
			}

			dist = sqrt(sum);

			if(min > dist) {
				min = dist;
				curr_point.label = i;
			}
		}
	}
}

vector<Point> averageLabelCentroids(vector<Point> points, vector<Cluster> clusters) {
	int total_points = points.size();

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < total_points; j++) {
			float sum = 0;
			int cluster_points = clusters[i].points.size();

			for (m = 0; m < total_points; m++) {
				sum += clusters[i].points[m].getCoordinate(j);
			}
			clusters[i]
		}
	}
}


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
