#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>

using namespace std;

int k;
float threshold;
int iterations;
int workers;
char* input;

class Point {
private:
	vector<float> coordinates;

public:

	Point(vector<float>& coordinates) {
		int size = coordinates.size();
		for(int i = 0; i < size; i++) {
			this->coordinates.push_back(coordinates[i]);
		}
	}

	float getDimensionsCount() {
		return coordinates.size();
	}

	float getCoordinate(int index) {
		return coordinates[index];
	}
};

class Centroid {
private:
	vector<Point> points;
	vector<float> coordinates;

public:
	Centroid(Point point) {
		int count = point.getDimensionsCount();
		for(int i = 0; i < count; i++) {
			coordinates.push_back(point.getCoordinate(i));
		}
		points.push_back(point);
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	float getCoordinate(int index) {
		return coordinates[index];
	}

	void setCoordinate(int index, float coordinate) {
		coordinates[index] = coordinate;
	}

	int getSize() {
		return points.size();
	}

	void addPoint(Point point)
	{
		points.push_back(point);
	}

	int getD() {
		return coordinates.size();
	}
};

vector<Centroid> randomCentroids(vector<Point> points, int k) {
	vector<Centroid> retval;
	int size = points.size();
	vector<int> cache;

	for (int i = 1; i < k + 1; i++) {
		LOOP:
		int x = rand() % (size - 1) + 1;
		if (find(cache.begin(), cache.end(), x) != cache.end())
			goto LOOP;
		Centroid centroid(points[x]);
		retval.push_back(centroid);
		// printf("%f\n", retval[i-1].getCoordinate(0));
	}

	return retval;
}


void findNearestCentroids(vector<Point> points, vector<Centroid> centroids) {
	int total_points = points.size();
	int total_coordinates = points[1].getDimensionsCount();
	printf("Total Points: %d\n", total_points);
	printf("Total Coordinates: %d\n", total_coordinates);

	for (int i = 0; i < total_points; i++) {
		Point curr_point = points[i];

		float sum_diff = 0;
		float min;

		for (int j = 0; j < total_coordinates; j++) {
			printf("Pow: %f %f\n", centroids[0].getCoordinate(j), curr_point.getCoordinate(j));
			sum_diff += pow(centroids[0].getCoordinate(j) - curr_point.getCoordinate(j), 2.0);
			printf("Sum: %f", sum_diff);
		}

		min = sqrt(sum_diff);

		for (int j = 1; j < k; j++) {
			float dist;
			float sum = 0;

			for(int m = 0; m < total_coordinates; m++) {
				sum += pow(centroids[j].getCoordinate(m) - curr_point.getCoordinate(m), 2.0);
			}

			dist = sqrt(sum);

			if(min > dist) {
				min = dist;
				centroids[j].addPoint(curr_point);
			}
		}
	}
}

vector<Centroid> averageLabeledCentroids(vector<Point> points, vector<Centroid> centroids) {
	int d = centroids[1].getD();
	int total_points = centroids[1].getSize();
	// printf("%d\n", total_points);
	vector<Centroid> retval;

	// printf("%f \n", centroids[0].getCoordinate(0));

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < d; j++) {
			float sum = 0;
			for (int m = 0; m < total_points; m++) {
				Point p = centroids[i].getPoint(m);
				// printf("%f + ", p.getCoordinate(j));
				sum = sum + p.getCoordinate(j);
			}
			// printf("= %f\n", sum);
			centroids[i].setCoordinate(j, sum / total_points);
		}
		retval.push_back(centroids[i]);
	}

	// printf("%f \n", centroids[0].getCoordinate(0));

	return retval;
}

bool converged(vector<Centroid> centroids, vector<Centroid> old_centroids) {
	int total = centroids.size();
	int d = centroids[0].getD();
	int diff;
	for (int i = 0; i < total; i++) {
		for (int j = 0; j < d; j++) {
			diff = abs(centroids[i].getCoordinate(j) - old_centroids[i].getCoordinate(j));
			if (diff > threshold)
				return false;
		}
	}
	return true;
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

		if (iterations > 0) {
			done = iters > iterations || converged(centroids, old_centroids);
		} else {
			done = converged(centroids, old_centroids);
		}
	}

	printf("Converged in %d iterations\n", iters);

	int d = centroids[0].getD();

	for (int i = 0; i < k; i++) {
		printf("Cluster %d center: ", i);
		for (int j = 0; j < d; j++) {
			printf("%f ", centroids[i].getCoordinate(j));
		}
		printf("\n");
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

	// Reading in input from a .txt file

	int rows;
	ifstream inFile;
	string line;
	string s;

	inFile.open("input.txt");

	inFile >> rows;

	// printf("Rows: %d\n", rows);

	vector<Point> dataset;
			vector<float> coordinates;
					int d;


	while (getline(inFile, line)) {
		istringstream iss(line);
		float x;
		iss >> x;
		d = 0;
		while(iss >> s)
        {
        	d++;
        	x = stof(s);
        	// printf("%f ", x);
            coordinates.push_back(x);
        }
        // printf("\n");
	}
	// printf("D: %d\n", d);
	Point point(coordinates);
    dataset.push_back(point);

    int i = d;
    while (i < d*rows) {
		vector<float> c;

    	for (int j = i; j < (i + d); j++) {
	    	c.push_back(dataset[0].getCoordinate(j));
	    }

	    Point point(c);
	    dataset.push_back(point);

	    i += 4;
    }

	Point x = dataset[0];
	float y = x.getCoordinate(0);
	// printf("Here\n");

	// for (int a = 0; a < rows; a++) {
	// 		printf("%f %f %f %f\n", dataset[a].getCoordinate(0),
	// 			dataset[a].getCoordinate(1), dataset[a].getCoordinate(2),
	// 			dataset[a].getCoordinate(3));
  //
	// }

	kmeans(dataset, k);
}
