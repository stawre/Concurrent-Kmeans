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
#include <time.h>
#include <pthread.h>

using namespace std;

int k;
float threshold;
int iterations;
int workers;
char* input;

pthread_barrier_t barrier;
pthread_mutex_t mutex;

class Point {
private:
	int id;
	int centroid_id = -1;
	vector<double> coordinates;

public:

	Point(int id, vector<double> coordinates) {
		this->id = id;
		int size = coordinates.size();
		for(int i = 0; i < size; i++) {
			this->coordinates.push_back(coordinates[i]);
		}
	}

	double getDimensionsCount() {
		return coordinates.size();
	}

	double getCoordinate(int index) {
		return coordinates[index];
	}

	int getId() {
		return id;
	}

	int getCentroid() {
		return centroid_id;
	}

	int setCentroid(int centroid_id) {
		this->centroid_id = centroid_id;
	}
};

class Centroid {
private:
	int id;
	vector<Point> points;
	vector<double> coordinates;

public:
	Centroid(int id, Point point) {
		this->id = id;
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

	double getCoordinate(int index) {
		return coordinates[index];
	}

	void setCoordinate(int index, double coordinate) {
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

	int getId() {
		return id;
	}

	bool erasePoint(int id)
	{
		int size = points.size();

		for(int i = 0; i < size; i++)
		{
			if(points[i].getId() == id)
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

};

vector<Centroid> randomCentroids(vector<Point>& points, int k) {
	vector<Centroid> retval;
	int size = points.size();
	vector<int> cache;
	int id = 0;
	for (int i = 1; i < k + 1; i++) {
		LOOP:
		int x = rand() % (size - 1) + 1;
		if (find(cache.begin(), cache.end(), x) != cache.end())
		goto LOOP;
		points[x].setCentroid(i - 1);
		Centroid centroid(id, points[x]);
		id++;
		retval.push_back(centroid);
		// printf("%f\n", retval[i-1].getCoordinate(0));
	}

	return retval;
}


void findNearestCentroids(vector<Point>& points, vector<Centroid>& centroids) {
	int total_points = points.size();
	int total_coordinates = points[1].getDimensionsCount();
	// printf("Total Points: %d\n", total_points);
	// printf("Total Coordinates: %d\n", total_coordinates);

	for (int i = 0; i < total_points; i++) {
		Point curr_point = points[i];

		double sum_diff = 0;
		double min;

		for (int j = 0; j < total_coordinates; j++) {
			// printf("Pow: %f %f\n", centroids[0].getCoordinate(j), curr_point.getCoordinate(j));
			sum_diff += pow(centroids[0].getCoordinate(j) - curr_point.getCoordinate(j), 2.0);
			// printf("Sum: %f\n", sum_diff);
		}

		min = sqrt(sum_diff);
		// printf("Min: %f\n", min);

		for (int j = 1; j < k; j++) {
			double dist;
			double sum = 0;

			for(int m = 0; m < total_coordinates; m++) {
				// printf("Pow: %f %f\n", centroids[j].getCoordinate(m), curr_point.getCoordinate(m));
				sum += pow(centroids[j].getCoordinate(m) - curr_point.getCoordinate(m), 2.0);
				// printf("Sum: %f\n", sum);
			}

			dist = sqrt(sum);
			// printf("Dist: %f\n", dist);

			if(min > dist) {
				min = dist;
				if (curr_point.getCentroid() != -1) {
					int temp = curr_point.getCentroid();
					// printf("Here\n");
					centroids[temp].erasePoint(curr_point.getId());
					// printf("Removed\n");
				}
				centroids[j].addPoint(curr_point);
				curr_point.setCentroid(j);
			} else {
				centroids[0].addPoint(curr_point);
				curr_point.setCentroid(0);
			}
			// printf("Final Min: %f\n", min);
		}
		// printf("Centroid size: %d\n", centroids[0].getSize());
	}
}

void averageLabeledCentroids(vector<Point>& points, vector<Centroid>& centroids) {
	int d = centroids[1].getD();
	//int total_points = centroids[1].getSize();
	//printf("Total Points: %d\n", total_points);
	// vector<Centroid> retval;

	// printf("%f \n", centroids[0].getCoordinate(0));

	for (int i = 0; i < k; i++) {
		int total_points = centroids[i].getSize();
		// printf("Total Points: %d\n", total_points);
		for (int j = 0; j < d; j++) {
			double sum = 0;
			for (int m = 0; m < total_points; m++) {
				Point p = centroids[i].getPoint(m);
				// printf("%f + ", p.getCoordinate(j));
				sum = sum + p.getCoordinate(j);
			}
			// printf("Coordinate: %f\n", sum / total_points);
			centroids[i].setCoordinate(j, sum / total_points);
		}
		// retval.push_back(centroids[i]);
	}

	// printf("%f \n", centroids[0].getCoordinate(0));

	// return retval;
}

bool converged(vector<Centroid> centroids, vector<Centroid> old_centroids) {
	int total = centroids.size();
	int d = centroids[0].getD();
	float diff;
	for (int i = 0; i < total; i++) {
		for (int j = 0; j < d; j++) {
			diff = abs(centroids[i].getCoordinate(j) - old_centroids[i].getCoordinate(j));
			// printf("Diff: %f\n", diff);
			if (diff > threshold)
			return false;
		}
	}
	return true;
}

void* kmeans(void* arg) {
	vector<Point>* dataset;
	dataset = (vector<Point> *) arg;
	vector<Centroid> centroids = randomCentroids(*dataset, k);

	int iters = 0;
	vector<Centroid> old_centroids;

	bool done = false;
	// for (int i = 0; i < k; i++) {
	//               printf("Cluster %d center: ", i);
	//               for (int j = 0; j < 4; j++) {
	//                       printf("%f ", centroids[i].getCoordinate(j));
	//               }
	//               printf("\n");
	//       }

	clock_t t;
	t = clock();

	while (!done) {
		pthread_barrier_wait(&barrier);
		old_centroids = centroids;
		iters++;

		pthread_mutex_lock(&mutex);
		findNearestCentroids(*dataset, centroids);
		pthread_mutex_unlock(&mutex);

		pthread_mutex_lock(&mutex);
		averageLabeledCentroids(*dataset, centroids);
		pthread_mutex_unlock(&mutex);

		pthread_barrier_wait(&barrier);

		if (iterations > 0) {
			done = iters > iterations || converged(centroids, old_centroids);
		} else {
			done = converged(centroids, old_centroids);
		}
	}

	for (int i = 0; i < dataset->size(); i++) {
		printf("Point %d label: %d\n", dataset[i]->getId(), dataset[i]->getCentroid());
	}

	t = clock() - t;
	double time_taken = ((double) t) / CLOCKS_PER_SEC;
	printf("Time taken by kmeans: %f seconds\n", time_taken);

	// printf("Centroid 0 size %d", centroids[0].getSize());
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
	vector<double> coordinates;
	int d;


	while (getline(inFile, line)) {
		istringstream iss(line);
		double x;
		iss >> x;
		d = 0;
		while(iss >> s)
		{
			d++;
			x = stof(s);
			// printf("%f ", x);
			coordinates.push_back(x);
		}
	}
	int id = 0;
	Point point(id, coordinates);
	id++;
	dataset.push_back(point);
	int i = d;
	while (i < d*rows) {
		vector<double> c;

		for (int j = i; j < (i + d); j++) {
			c.push_back(dataset[0].getCoordinate(j));
		}

		Point point(id, c);
		id++;
		dataset.push_back(point);

		i += 4;
	}

	// printf("Dataset size: %d\n", dataset.size());

	pthread_t threads[workers];
	pthread_barrier_init(&barrier, NULL, workers);
	pthread_mutex_init(&mutex, NULL);

	for (i = 0; i < workers; i++)
	{
		pthread_create(&threads[i], NULL, kmeans, &dataset);
	}

	for (i = 0; i < workers; i++)
	{
		pthread_join(threads[i], NULL);
	}

	// kmeans(dataset, k);
}
