#include "Network.h"
#include "Common.h"

#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

double Network::gauss(int dist, double var) {
	return exp(-0.5*pow((dist/var), 2));
}

int Network::calcXCoord(int width, int index) {
	return index % width;
}

int Network::calcYCoord(int width, int index) {
	return (index / width);
}

int Network::calcIndex(int width, int x, int y) {
	return (width*y + x);
}

double Network::calcDistance(int width, int pre, int post) {
	int pre_x = calcXCoord(width, pre);
	int pre_y = calcYCoord(width, pre);
	int post_x = calcXCoord(width, post);
	int post_y = calcYCoord(width, post);

	return sqrt(pow(pre_x - post_x, 2.0) + pow(pre_y - post_y, 2.0));
}

Network::Network(struct network_constants_s *network_constants,
		struct layer_constants_s *layer_constants) {

	/* Copy network constants */
	nc.width = network_constants->width;
	nc.num_layers = network_constants->num_layers;
	nc.radius = network_constants->radius;

	/* Allocate and copy layer constants */
	lc = new struct layer_constants_s[nc.num_layers-1];
	memcpy(lc, layer_constants, (nc.num_layers-1)*sizeof(struct layer_constants_s));

	/* Initialise stat vars */
	num_neurons = nc.num_layers*nc.width*nc.width;
	num_synapses = 0;

	/* Allocate neurons */
	neurons = new double[num_neurons];
	
	/* Initialise neurons */
	for (int i=0; i<num_neurons; i++)
		neurons[i] = 0.0;

	/* Create synapse deques */
	synapses = new deque<struct synapse_s*>[nc.num_layers-1];

	/* Generate synapses */
	srand(time(NULL));
	for (int l=0; l<nc.num_layers-1; l++) {
		for (int x=0; x<nc.width*nc.width; x++) {
			for (int y=0; y<nc.width*nc.width; y++) {
				double dist = calcDistance(nc.width, x, y);
				double pr = gauss(dist, nc.radius*2);
				double dr = (double)rand() / RAND_MAX;

				if (dr < pr) {
					/* Create synapse */
					struct synapse_s *syn = new struct synapse_s;
					syn->pre = x;
					syn->post = y;
					syn->weight = ((double)rand() / RAND_MAX);

					/* Add synapse */
					synapses[l].push_back(syn);

					num_synapses++;
				}
			}
		}
	}

	/* DEBUG */
	/*
	int n_x = nc.width/2;
	int n_y = n_x;
	int neuron = nc.width*n_y + n_x;
	cout << "Neuron (" << n_x << ", " << n_y
		<< ") has connections: ";
	for (int s=0; s<synapses[0].size(); s++) {
		if (synapses[0][s]->post == neuron)
			cout << calcXCoord(nc.width, synapses[0][s]->pre)
				<< " " << calcYCoord(nc.width, synapses[0][s]->pre)
				<< endl;
	}
	cout << endl;
	*/
}

void Network::applyNoise(int layer) {
	int width = nc.width;

	/* Generate and apply noise to top layer neurons */
	//for (int i=0; i<width*width; i++)
	//	neurons[i] = ((double)rand() / RAND_MAX)*2-1;

	/* Apply noise with low(er) resolution */
	int n = width/4;
	for (int n1=0; n1<n; n1++) {
		for (int n2=0; n2<n; n2++) {
			double r = ((double)rand() / RAND_MAX)*2-1;
			for (int w1=0; w1<4; w1++) {
				for (int w2=0; w2<4; w2++) {
					int i = calcIndex(width, n2*4 + w1, n1*4 + w2);
					neurons[i] = r;
				}
			}
		}
	}

	/* Zero all other neurons */
	for (int n=width*width; n<num_neurons; n++)
		neurons[n] = 0.0;

	for (int l=0; l<nc.num_layers-1; l++) {
		/* DEBUG */
		cout << "Computing layer " << l+1 << " energies... ";
		cout.flush();

		/* Propagate energies down a layer */
		for (int s=0; s<synapses[l].size(); s++) {
			int pre = synapses[l][s]->pre;
			int post = synapses[l][s]->post;
			double weight = synapses[l][s]->weight;

			/*
			Increment the energy by the energy of the
			presynaptic neuron multiplied by the synapse weight
			*/
			double energy_d = neurons[l*width*width + pre] * weight;
			neurons[(l+1)*width*width + post] += energy_d;

			/* DEBUG */
			//cout << "Added energy delta " << energy_d <<
			//	" from neuron " << pre << " to neuron " << post <<
			//	" in layer " << l+1 << endl;
		}

		/* Apply layer constants to summation term */
		for (int n=0; n<width*width; n++) {
			neurons[(l+1)*width*width + n] *= lc[l].Rb;
			neurons[(l+1)*width*width + n] += lc[l].Ra;
		}

		cout << "done." << endl;
		cout.flush();
	}

	/* Adjust synapses for 'layer' */
	cout << "Adjusting synapses... ";
	cout.flush();
	for (int s=0; s<synapses[layer].size(); s++) {
		/* Is weight saturated? */
		/* TODO move saturated weights into new deque (efficient) */
		if (synapses[layer][s]->weight == 1.0 ||
			synapses[layer][s]->weight == -1.0)
			continue;

		double b = lc[layer].b;
		double c = lc[layer].c;
		int pre = synapses[layer][s]->pre;
		int post = synapses[layer][s]->post;
		double pre_e = neurons[pre];
		double post_e = neurons[post];

		double weight_d = b + c*pre_e*post_e;
		double weight_new = synapses[layer][s]->weight + weight_d;

		if (weight_new > 1.0)
			weight_new = 1.0;
		else if (weight_new < -1.0)
			weight_new = -1.0;

		/* Change weight */
		synapses[layer][s]->weight = weight_new;

		/* DEBUG */
		//cout << "Weight delta: " << weight_d << ", " << 
		//	"new weight: " << weight_new << endl;
	}
	cout << "done." << endl;
	cout.flush();
}

void Network::writeSynapseData(int l, int x, int y, ofstream &outStream) {
	int index = calcIndex(nc.width, x, y);
	int pre_i, pre_x, pre_y;
	double weight;

	for (int s=0; s<synapses[l-1].size(); s++) {
		/* If this synapse has (x, y) as its post-syn neuron, print */
		if (synapses[l-1][s]->post == index) {
			pre_i = synapses[l-1][s]->pre;
			pre_x = calcXCoord(nc.width, pre_i);
			pre_y = calcYCoord(nc.width, pre_i);
			weight = synapses[l-1][s]->weight;

			outStream << pre_x << " " << pre_y << " " << weight << endl;
		}
	}

	cout.flush();
}

Network::~Network() {
	/* Delete layer constants */
	delete[] lc;

	/* Delete neurons */
	delete[] neurons;

	/* Delete synapses */
	for (int l=0; l<nc.num_layers-1; l++) {
		for (int i=0; i<synapses[l].size(); i++)
			delete synapses[l][i]; /* Delete each struct synapse_s */
	}
	delete[] synapses; /* Delete the array of deques */
}