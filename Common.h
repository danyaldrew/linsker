#ifndef H_COMMON
#define H_COMMON

struct synapse_s {
	int pre, post;
	double weight;
};

struct network_constants_s {
	int width, num_layers;
	double radius;
};

struct layer_constants_s {
	double b, c, Ra, Rb;
};

#endif