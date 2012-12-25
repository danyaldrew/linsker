#include <iostream>
#include <fstream>

#include "Common.h"
#include "Network.h"

using namespace std;

int main(int argc, char **argv) {
	struct network_constants_s nc = { 128, 3, 3.0 };
	struct layer_constants_s lc[] = { { 0.8, -0.3, 0.0, 1.0 },
		{ 0.01, -0.9, 0.0, 1.0 } };

	cout << "Generating network... "; cout.flush();
	Network *n = new Network(&nc, lc);
	cout << "done." << endl; cout.flush();

	/* Develop first layer */
	for (int i=0; i<5; i++) {
		ofstream dataFile;
		char filename[16];
		snprintf(filename, 16, "data.b.%i", i);

		dataFile.open(filename);

		n->applyNoise(0);
		n->writeSynapseData(1, 64, 64, dataFile);

		dataFile.close();
	}

	/* Develop second layer */
	for (int i=0; i<100; i++) {
		ofstream dataFile;
		char filename[16];
		snprintf(filename, 16, "data.c.%i", i);

		dataFile.open(filename);

		n->applyNoise(1);
		n->writeSynapseData(2, 64, 64, dataFile);

		dataFile.close();
	}

	cout << "Cleaning up... "; cout.flush();
	delete n;
	cout << "done." << endl; cout.flush();

	return 0;
}