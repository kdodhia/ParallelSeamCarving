#ifndef COMPUTE_ENERGY_ENTROPY_H
#define COMPUTE_ENERGY_ENTROPY_H

#define WINDOW 4


float compute_entropy(int *input, int size);

void get_entropy_map(int *image_gray, float *entropy_map, int width, int height);

#endif 