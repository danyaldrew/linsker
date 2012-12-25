#include <deque>
#include <fstream>

#include "Common.h"

using namespace std;

class Network {
private:
	/* Constants */
	struct network_constants_s nc;
	struct layer_constants_s *lc;

	/* Network Data */
	double *neurons;
	std::deque<struct synapse_s*> *synapses;

	/* Stat variables */
	int num_neurons, num_synapses;

	static double gauss(int dist, double var);
	static int calcXCoord(int width, int index);
	static int calcYCoord(int width, int index);
	static int calcIndex(int width, int x, int y);
	static double calcDistance(int width, int pre, int post);

public:
	Network(struct network_constants_s *network_constants,
			struct layer_constants_s *layer_constants);

	void applyNoise(int layer);
	void writeSynapseData(int l, int x, int y, ofstream &outStream);
	
	~Network();
};