#include <math.h>
#include <set>
#include <iterator>
#include <algorithm>
#include "compute_energy_entropy.h"

#include <opencv/cv.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <sys/time.h>

using namespace std;
using namespace cv;

float compute_entropy(int *input, int size){

	std::set<int> symset (input, input+size);
	int set_size = symset.size();
	float prob[set_size];
	
	std::set<int>::iterator iter, start, end;
	start = symset.begin();
	end = symset.end();
	int i = 0;

	for (iter = start; iter != end; ++iter){
		int n = std::count(input, input+size, *iter);
		prob[i] = n/(1.0 * size);
		i++;
	}

	float entropy = 0;
	for (int j = 0; j < set_size; j++){
		int p = prob[j];
		int ent = - (p * logb(p));
		entropy += ent;
	}

	return entropy;
} 

//
void get_entropy_map(Mat image_gray, float *entropy_map, int height, int width){

	int w_size = (2 * WINDOW) + 1;

	for (int row = 0; row < height; row++){

		for (int col = 0; col < width; col++){

			int region[w_size * w_size];
			int region_size = 0;
			int start_x = MAX(row - WINDOW, 0);
			int start_y = MAX(col - WINDOW, 0);
			int end_x = MIN(row + WINDOW, height);
			int end_y = MIN(col + WINDOW, width);

			for (int i = start_y; i < end_y + 1; i++){

				for (int j = start_x; j < end_x + 1; j++){

					//int img_index = i * height + j;
					region[(i-start_y)*(w_size) + j] = image_gray.at<uchar>(i, j);
					region_size++;
				}
			}

			float entropy = compute_entropy(region, region_size);
			entropy_map[row * width + col] = entropy;
		}
	}
}
