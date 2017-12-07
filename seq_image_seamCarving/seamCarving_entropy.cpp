#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>

#include <pthread.h>

#include <opencv/cv.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <sys/time.h>

#include "compute_energy_entropy.h"

using namespace std;
using namespace cv;


typedef unsigned long long timestamp_t;
void calculate_energy(Mat greyImg, int *energy, int rows, int cols);
extern void get_entropy_map(Mat image_gray, float *entropy_map, 
                            int width, int height);
Mat reduce_image(Mat image,  int v, int h);
void calculate_ACM(int *energy, int rows, int cols);
void find_seam(int *energy, int *seam, int rows, int cols);
Mat remove_seam(Mat greyImg, Mat &image, int *seam, int rows, int cols);

timestamp_t get_timestamp () {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int main(int argc, char* argv[]) {
    cout << "Program started" << endl;
    timestamp_t t0 = get_timestamp();

    int ver = 0;
    int hor = 0;

    string inFile = "image.jpg";
    if(argc == 4) {
        string inFile = argv[1];
        ver = atoi(argv[2]);
        hor = atoi(argv[3]);
    }
    else {
        return -1;
    }

    Mat image = imread(inFile, CV_LOAD_IMAGE_COLOR); 
    if(!image.data){
        cerr << "Unable to open input file." << endl;
        return -1;
    }

    cout << "Number of cols before: " << image.cols << endl;

    string::size_type pAt = inFile.find_last_of('.');   // Find extension point
    const string outFile = inFile.substr(0, pAt) + "-result.jpg";
    Mat outImage = reduce_image(image, ver, hor);
    imwrite(outFile, outImage);

    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "Total runtime: " << secs << "secs" << endl;
    cout << "Number of cols after: " << outImage.cols << endl;

    return 0;
}


void calculate_energy(Mat greyImg, float *energy, int rows, int cols){
    cout << "generating energy matrix" << endl;
    timestamp_t t0 = get_timestamp();
    
    for(int col = 0; col < cols; col++){
        for(int row = 1; row < rows; row++){
            float energy_val = 0;
            
            if(col > 0 && col + 1 < cols)
                energy_val += abs(greyImg.at<uchar>(row,col + 1) - greyImg.at<uchar>(row,col - 1));
            else if(col > 0)
                energy_val += 2 * abs(greyImg.at<uchar>(row,col) - greyImg.at<uchar>(row,col - 1));
            else
                energy_val += 2 * abs(greyImg.at<uchar>(row,col + 1) - greyImg.at<uchar>(row,col));
            
            if(row > 0 && row + 1 < rows)
                energy_val += abs(greyImg.at<uchar>(row + 1,col) - greyImg.at<uchar>(row - 1,col));
            else if(row > 0)
                energy_val += 2 * abs(greyImg.at<uchar>(row,col) - greyImg.at<uchar>(row - 1,col));
            else
                energy_val += 2 * abs(greyImg.at<uchar>(row + 1,col) - greyImg.at<uchar>(row,col));
            
            energy[row * cols + col] = energy_val;
        }
    }
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "Energy matrix generation time: " << secs << "secs" << endl;
}

void find_seam(float *energy, int *seam, int rows, int cols)
{
    cout << "finding seam" << endl;
    timestamp_t t0 = get_timestamp();
    
    float  min_val = -1;
    int min_col = -1;
    int inf = 100000;

    for (int col = 0; col < cols-1; col++) {
        float cur = energy[(rows-1) * cols + col];
        if (cur < min_val || min_col == -1 ) {
            min_val = cur;
            min_col = col;
        }
    }

    int cur_col = min_col;
    int row;
    for (row = rows-1; row > 0; row--) {
        seam[row] = cur_col;
        float  upleft, upright;
        if (cur_col > 1) {
            upleft = energy[(row-1) * cols + (cur_col-1)];
        } else {
            upleft = inf;
        }

        if (cur_col < cols-2) {
            upright = energy[(row-1) * cols + (cur_col+1)];
        } else {
            upright = inf;
        }

        float up = energy[(row-1)*cols + cur_col];

        min_val = min(upleft, min(up, upright));

        if (min_val == upright) {
            cur_col++;
        } else if (min_val == upleft) {
            cur_col--;
        }
    }
    seam[row] = cur_col;

    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "Seam generation time: " << secs << "secs" << endl;
}

/*
 * Remove seam from the image.
 */
Mat remove_seam(Mat greyImg, Mat &image, int *seam, int rows, int cols)
{
    cout << "Generating ACM" << endl;
    timestamp_t t0 = get_timestamp();

    Mat reducedImage(rows,cols-1,CV_8UC3, Scalar(0, 0, 0));
    Mat reducedImageGrey(rows,cols-1,CV_8UC1);
    for (int row = 0; row < rows; row++) {
        int col_to_remove = seam[row];
        for (int col = 0; col < cols; col++) {
            if (col > col_to_remove) {
                reducedImage.at<Vec3b>(row, col-1) = image.at<Vec3b>(row, col);
                reducedImageGrey.at<uchar>(row, col-1) = greyImg.at<uchar>(row, col);
            } else if (col < col_to_remove) {
                reducedImage.at<Vec3b>(row, col) = image.at<Vec3b>(row, col);
                reducedImageGrey.at<uchar>(row, col) = greyImg.at<uchar>(row, col);
            }
        }
    }
    image = reducedImage;
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "ACM generation time: " << secs << "secs" << endl;

    return reducedImageGrey;
}

void calculate_ACM(float *energy, int rows, int cols) {
        
    cout << "Generating ACM" << endl;
    timestamp_t t0 = get_timestamp();

    for (int row = 2; row < rows; row++) {
        for (int col = 1; col < cols; col++) {
            float up = energy[(row-1) * cols + col];
            if (col == 1) {
                float upright = energy[(row-1) * cols + (col+1)];
                energy[row * cols + col] = energy[row * cols + col] + min(up, upright);
            } else if (col == cols - 2){
                float upleft = energy[(row-1) * cols + col-1];
                energy[row * cols + col] = energy[row * cols + col] + min(up, upleft);
            } else {
                float upright = energy[(row-1) * cols + (col+1)];
                float upleft = energy[(row-1) * cols + (col-1)];
                energy[row * cols + col] = energy[row * cols + col] + min(up, min(upleft, upright));
            }
        }
    }
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "ACM generation time: " << secs << "secs" << endl;
}



Mat reduce_image(Mat image,  int v, int h) {

    Mat grayImg1;
    cvtColor(image, grayImg1, CV_RGB2GRAY);

    int min = 0, diff = 0;
    Mat reducedGrayImg1, reducedImg;
    grayImg1.copyTo(reducedGrayImg1);
    
    image.copyTo(reducedImg);

    int cols = image.cols;
    int rows = image.rows;

    float *energy = (float *)calloc(cols * rows, sizeof(float));
    int *seam = (int *)calloc(rows, sizeof(int));

    for(int i = 0; i < v; ++i) {
        /*
         * Remove one vertical seam from img. The algorithm:
        1) get energy matrix.
        2) generate accumulated energy matrix
        3) find seam
        4) remove this seam from the image
         */

        //calculate_energy(reducedGrayImg1, energy, rows, cols);
        get_entropy_map(reducedGrayImg1, energy, rows, cols);

        calculate_ACM(energy, rows, cols);

        find_seam(energy, seam, rows, cols);

        reducedGrayImg1 = remove_seam(reducedGrayImg1, reducedImg, seam, rows, cols);

        cols--;
    }
    return reducedImg;
}