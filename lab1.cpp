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
double threshold;
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

vector<Point> dataset;
vector<Centroid> centroids;
vector<vector<Point>> partitions;

typedef struct thread_data {
		int i;
		vector<Point> partition;
} data_t;

void randomCentroids() {
	vector<Centroid> retval;
	int size = dataset.size();
	vector<int> cache;
	int id = 0;
	for (int i = 1; i < k + 1; i++) {
		LOOP:
		// printf("Here 0\n");
		// printf("Size: %d\n", size);
		int x = rand() % size + 1;
		// printf("Rand: %d\n", x);
		// printf("Here 1\n");

		if (find(cache.begin(), cache.end(), x) != cache.end())
			goto LOOP;
		cache.push_back(x);
		dataset[x].setCentroid(i - 1);
		Centroid centroid(id, dataset[x]);
		id++;
		centroids.push_back(centroid);
		// printf("%f\n", retval[i-1].getCoordinate(0));
	}

	// return retval;
}


void* findNearestCentroids(void *arg) {
	int a;
	vector<Point> dataset ;
	data_t* my_data;

	printf("Hello 0\n");

	my_data = (data_t *) arg;
	a = my_data->i;
	printf("i \n");
	dataset = my_data->partition;

	printf("Hey \n");

	int total_points = dataset.size();
	int total_coordinates = dataset[1].getDimensionsCount();
	// printf("Total Points: %d\n", total_points);
	// printf("Total Coordinates: %d\n", total_coordinates);

	printf("Hello 1\n");

	for (int i = 0; i < total_points; i++) {
		double sum_diff = 0;
		double min;
		int index = 0;

		printf("Hello 2\n");

		for (int j = 0; j < total_coordinates; j++) {
			pthread_mutex_lock(&mutex);
			sum_diff += pow(centroids[0].getCoordinate(j) - dataset[i].getCoordinate(j), 2.0);
			pthread_mutex_unlock(&mutex);
		}

		printf("Hello 3\n");

		min = sqrt(sum_diff);

		for (int j = 1; j < k; j++) {
			double dist;
			sum_diff = 0;

			for(int m = 0; m < total_coordinates; m++) {
				pthread_mutex_lock(&mutex);
				sum_diff += pow(centroids[j].getCoordinate(m) - dataset[i].getCoordinate(m), 2.0);
				pthread_mutex_unlock(&mutex);
			}

			dist = sqrt(sum_diff);

			if(min > dist) {
				min = dist;
				index = j;				
			}
		}
		// printf("Index: %d\n", index);
		// printf("Centroid size: %d\n", centroids[0].getSize());
		if (index != dataset[i].getCentroid()) {
			if (dataset[i].getCentroid() != -1) {
				pthread_mutex_lock(&mutex);
				int c_id = dataset[i].getCentroid();
				// printf("Here\n");
				centroids[c_id].erasePoint(dataset[i].getId());
				pthread_mutex_unlock(&mutex);
				// printf("Removed\n");
			}
	
			pthread_mutex_lock(&mutex);
			centroids[index].addPoint(dataset[i]);
			dataset[i].setCentroid(index);
			pthread_mutex_unlock(&mutex);
		}
	}
	pthread_mutex_lock(&mutex);
	partitions[a] = dataset;
	pthread_mutex_unlock(&mutex);
}

void averageLabeledCentroids() {
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

bool converged(vector<Centroid> old_centroids) {
	int total = centroids.size();
	int d = centroids[0].getD();
	double diff;
	for (int i = 0; i < total; i++) {
		for (int j = 0; j < d; j++) {
			diff = abs(centroids[i].getCoordinate(j) - old_centroids[i].getCoordinate(j));
			printf("Diff: %f\n", diff);
			if (diff > threshold)
				return false;
		}
	}
	return true;
}

void parallel_solution() {
	int total_points = dataset.size();
	int per_part = total_points / workers;

	int counter = 0;


	for (int i = 0; i < workers; i++) {
		vector<Point> p;
		for (int j = 0; j < per_part; j++) {
			p.push_back(dataset[counter++]);
		}
		partitions.push_back(p);
	}
	
	int rem = total_points % workers;

	if (rem != 0 && counter < total_points) {
		for (int i = 0; i < rem; i++) {
			partitions[i].push_back(dataset[counter++]);
		}
	}

	pthread_t threads[workers];
	data_t t_data[workers];
	
	for (int i = 0; i < workers; i++) {
		t_data[i].i = i;
		t_data[i].partition = partitions[i];
	}

	for (int i = 0; i < workers; i++) {
		pthread_create(&threads[i], NULL, findNearestCentroids, &t_data[i]);
	}

	for (int i = 0; i < workers; i++)
	{
		pthread_join(threads[i], NULL);
	}
}

void kmeans() {
	// vector<Point>* dataset;
	// dataset = (vector<Point> *) arg;
	randomCentroids();

	int iters = 0;
	vector<Centroid> old_centroids;

	bool done = false;

	clock_t t;
	t = clock();	

	while (!done) {
		// pthread_barrier_wait(&barrier);
		old_centroids = centroids;
		iters++;

		// pthread_mutex_lock(&mutex);
		printf("Here 0\n");
		parallel_solution();
		// pthread_mutex_unlock(&mutex);
		printf("Here 1\n");

		// pthread_mutex_lock(&mutex);
		averageLabeledCentroids();
		// pthread_mutex_unlock(&mutex);

		pthread_barrier_wait(&barrier);
		if (iterations > 0) {
			done = iters > iterations || converged(old_centroids);
		} else {
			done = converged(old_centroids);
		}
	}

	vector<Point> data_copy = dataset;

	// for (int i = 0; i < data_copy.size(); i++) {
	// 	printf("Point %d label: %d\n", data_copy[i].getId(), data_copy[i].getCentroid());
	// }

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

	inFile.open(input);

	inFile >> rows;

	// printf("Rows: %d\n", rows);

	// vector<Point> dataset;
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
	int i = 32;
	while (i < 2097152) {
		vector<double> c;

		for (int j = i; j < (i + 32); j++) {
			c.push_back(dataset[0].getCoordinate(j));
		}

		Point point(id, c);
		id++;
		dataset.push_back(point);

		i += 32;
	}

	// for (int i = 0; i < dataset.size(); i++) {
	// 	printf("Point %d label: %d\n", dataset[i].getId(), dataset[i].getCentroid());
	// }

	// for (int i = 0; i < dataset.size(); i++) {
	// 	printf("Point %d: ", i);
	// 	for (int j = 0; j < 9; j++) {
	// 		printf("%f ", dataset[i].getCoordinate(j));
	// 	}
	// 	printf("\n");
	// }

	// printf("Dataset size: %d\n", dataset.size());

	pthread_barrier_init(&barrier, NULL, workers);
	pthread_mutex_init(&mutex, NULL);

	kmeans();
}
