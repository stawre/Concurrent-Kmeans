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
// pthread_spinlock_t spin;
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

	void setPoint(int index, Point point)
	{
		if (index < points.size()) {
			printf("In Set Point\n");
			points[index] = point;
		} else {
			points.push_back(point);
		}
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
vector<Centroid> old_centroids;
vector<vector<Centroid>> private_centroids;
vector<vector<Point>> partitions;

int sizes[16][16] = {0};

void populate_privates() {
	// int rem = k % workers;
	// int quotient = k / workers;
	// int counter = 0;

	vector<Centroid> centroid;

	for (int i = 0; i < workers; i++) {
		private_centroids.push_back(centroid);
	}

	printf("Hello -1\n");
	
	for (int i = 0; i < workers; i++) {
		for (int j = 0; j < k; j++) {
			printf("Hello 0\n");
			printf("%d %d\n", j, k);
			private_centroids[i].push_back(centroids[j]);
			// counter++;
		}
	}
	printf("Hello 1\n");		
}	

void integrate_centroids() {
	int rem = k % workers;
        int quotient = k / workers;
	int counter = 0;

	for(int i = 0; i < workers; i++) {
		for(int j = 0; j < quotient; j++) {
			centroids[counter++] = private_centroids[i][j];
		} 
	}
}


typedef struct thread_data {
	int i;
	vector<Point> partition;
} data_t;

void randomCentroids() {
	// vector<Centroid> retval;
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

	old_centroids = centroids;

	// return retval;
}


void findNearestCentroids(vector<Point>& my_dataset, vector<Centroid> old_centroids2, int a) {
	int total_points = my_dataset.size();
	int total_coordinates = my_dataset[1].getDimensionsCount();
	// printf("Total Points: %d\n", total_points);
	// printf("Total Coordinates: %d\n", total_coordinates);

	// printf("Hello 1\n");

	for (int i = 0; i < total_points; i++) {
		double sum_diff = 0;
		double min;
		int index = 0;

		// printf("Hello 2\n");

		for (int j = 0; j < total_coordinates; j++) {
			// pthread_mutex_lock(&mutex);
			sum_diff += pow(old_centroids2[0].getCoordinate(j) - my_dataset[i].getCoordinate(j), 2.0);
			// pthread_mutex_unlock(&mutex);
		}

		// printf("Hello 3\n");

		min = sqrt(sum_diff);

		for (int j = 1; j < k; j++) {
			double dist;
			sum_diff = 0;
			// printf("Hello 4\n");
			for(int m = 0; m < total_coordinates; m++) {
				// printf("%f\n", old_centroids[j].getCoordinate(m));
				// pthread_mutex_lock(&mutex);
				sum_diff += pow(old_centroids2[j].getCoordinate(m) - my_dataset[i].getCoordinate(m), 2.0);
				// pthread_mutex_unlock(&mutex);
			}

			dist = sqrt(sum_diff);

			if(min > dist) {
				min = dist;
				index = j;
			}
		}
		// printf("Hello 5\n");
		// printf("Index: %d\n", index);
		// printf("Centroid size: %d\n", centroids[0].getSize());
		if (index != my_dataset[i].getCentroid()) {
			if (my_dataset[i].getCentroid() != -1) {
				// pthread_mutex_lock(&mutex);
				int c_id = my_dataset[i].getCentroid();
				// printf("Here\n");
				centroids[c_id].erasePoint(my_dataset[i].getId());
				// pthread_mutex_unlock(&mutex);
				// printf("Removed\n");
			}

			// pthread_mutex_lock(&mutex);
			old_centroids2[index].addPoint(my_dataset[i]);
			sizes[a][index]++;
			// pthread_mutex_unlock(&mutex);
			my_dataset[i].setCentroid(index);
		}
	}
	printf("Hey\n");
	// pthread_spin_lock(&spin);
	printf("Hey there\n");
	private_centroids[a] = old_centroids2;
	// pthread_spin_unlock(&spin);
	// printf("Hey there\n");
	// pthread_mutex_lock(&mutex);
	partitions[a] = my_dataset;
	// pthread_mutex_unlock(&mutex);
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
		// printf("Average 1\n");
		for (int j = 0; j < d; j++) {
			double sum = 0;
			// printf("Average 2\n");
			for (int m = 0; m < total_points; m++) {
				Point p = centroids[i].getPoint(m);
				// printf("%f + ", p.getCoordinate(j));
				sum = sum + p.getCoordinate(j);
			}
			// printf("Average 3\n");
			// printf("Coordinate: %f\n", sum / total_points);
			centroids[i].setCoordinate(j, sum / total_points);
		}
		// retval.push_back(centroids[i]);
		// printf("Average 4\n");
	}

	// printf("%f \n", centroids[0].getCoordinate(0));

	// return retval;
}

bool converged() {
	int total = centroids.size();
	int d = centroids[0].getD();
	double diff;
	for (int i = 0; i < total; i++) {
		// printf("Converged 1\n");
		for (int j = 0; j < d; j++) {
			// printf("%f\n", old_centroids[i].getCoordinate(j));
			diff = abs(centroids[i].getCoordinate(j) - old_centroids[i].getCoordinate(j));
			// printf("Diff: %f\n", diff);
			if (diff > threshold)
			return false;
		}
		// printf("Converged 2\n");
	}
	return true;
}

void* kmeans(void *arg) {
	// vector<Point>* dataset;
	// dataset = (vector<Point> *) arg;
	int a;
	vector<Point> t_dataset ;
	data_t* my_data;

	my_data = (data_t *) arg;
	a = my_data->i;
	t_dataset = my_data->partition;

	// randomCentroids();

	int iters = 0;
	bool done = false;

	clock_t t;
	t = clock();

	while (!done) {
		pthread_barrier_wait(&barrier);
		if (a == 0) {
			old_centroids = centroids;
		 	iters++;
		}
		// old_centroids = private_centroids[a];
		// pthread_barrier_wait(&barrier);
		// pthread_mutex_lock(&mutex);
		printf("Here 0\n");
		findNearestCentroids(t_dataset, old_centroids, a);
		// pthread_mutex_unlock(&mutex);
		printf("Here 1\n");

		// pthread_barrier_wait(&barrier);
		// pthread_mutex_lock(&mutex);
		// integrate_centroids();
		// pthread_mutex_unlock(&mutex);
		
		pthread_barrier_wait(&barrier);
		pthread_mutex_lock(&mutex);

		for (int i = 0; i < k; i++) {
			// int index = indexes[a];
			for (int j = 0; j < sizes[a][i]; j++) {
				centroids[i].addPoint(private_centroids[a][i].getPoint(j));
			}
		}
		pthread_mutex_unlock(&mutex);

		pthread_barrier_wait(&barrier);

		if (a == 0)
			averageLabeledCentroids();
		// pthread_mutex_unlock(&mutex);

		printf("Here 2\n");

		pthread_barrier_wait(&barrier);
		// printf("Here 3\n");
		// pthread_mutex_lock(&mutex);
		if (a == 0) {
			if (iterations > 0) {
				// printf("Here 4\n");
				done = iters > iterations || converged();
			} else {
				done = converged();
			}
		}
		// pthread_mutex_unlock(&mutex);
		// pthread_barrier_wait(&barrier);
	}

	// vector<Point> data_copy = dataset;

	// for (int i = 0; i < data_copy.size(); i++) {
	// 	printf("Point %d label: %d\n", data_copy[i].getId(), data_copy[i].getCentroid());
	// }

	if (a == 0) {
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
		pthread_create(&threads[i], NULL, kmeans, &t_data[i]);
	}

	for (int i = 0; i < workers; i++)
	{
		pthread_join(threads[i], NULL);
	}

	counter = 0;

	for (int i = 0; i < workers; i++) {
		for (int j = 0; j < per_part; j++) {
			dataset[counter++] = partitions[i][j];
		}
	}

	if (rem != 0 && counter < total_points) {
		for (int i = 0; i < rem; i++) {
			dataset[counter++] = partitions[i][per_part];
		}
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
	//pthread_spin_init(&spin, 0);
	pthread_mutex_init(&mutex, NULL);

	randomCentroids();
	printf("Hey 1\n");
	populate_privates(); // seperate global 2D vector to store private updates
	printf("Hey 2\n");
	parallel_solution();
}
