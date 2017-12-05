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

using namespace std;
using namespace cv;

int ver = 1;
int hor = 0;
Mat image;
Mat outImage;
typedef unsigned long long timestamp_t;

void usage();
Mat reduce();
Mat remove_seam(Mat image, int seam[]);
Mat remove_seam_gray(Mat GrayImage, int seam[]);
int *find_seam(Mat grayImg1);
Mat reduce_vertical(Mat &grayImg1, Mat img);
Mat reduce_horizontal(Mat &grayImg1, Mat img);
Mat reduce_frame(Mat frame1, int v, int h);
Mat get_energy_matrix(Mat grayImg);
Mat generate_ACM(Mat energy);
double min_helper(double x, double y, double z);

timestamp_t get_timestamp () {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int main(int argc, char* argv[]) {
    cout << "Program started" << endl;
    timestamp_t t0 = get_timestamp();
    string inFile = "image.jpg";
    if(argc == 4) {
        string inFile = argv[1];
        ver = atoi(argv[2]);
        hor = atoi(argv[3]);
    }
    else {
        usage();
        return -1;
    }
    image = imread(inFile, CV_LOAD_IMAGE_COLOR); 
    if(!image.data){
        cerr << "Unable to open input file." << endl;
        return -1;
    }
    cout << "Number of cols before: " << image.cols << endl;

    string::size_type pAt = inFile.find_last_of('.');   // Find extension point
    const string outFile = inFile.substr(0, pAt) + "-result.jpg";
    outImage = reduce();
    imwrite(outFile, outImage);
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "Total runtime: " << secs << "secs" << endl;
    cout << "Number of cols after: " << outImage.cols << endl;
    return 0;
}


void usage() 
{
    cout << "Usage: heuristic5 <filename> <vertical cuts> <horizontal cuts>" << endl;
}

Mat reduce() {
    Mat frame = image;
    outImage = reduce_frame(frame, ver, hor);
    return outImage;
}



/*
 * Remove seam from the image.
 */
Mat remove_seam(Mat image, int seam[])
{
    int rows = image.rows;
    int cols = image.cols;
    Mat reducedImage(rows,cols-1,CV_8UC3, Scalar(0, 0, 0));
    for (int row = 0; row < rows; row++) {
        int col_to_remove = seam[row];
        for (int col = 0; col < cols; col++) {
            if (col > col_to_remove) {
                reducedImage.at<Vec3b>(row, col-1) = image.at<Vec3b>(row, col);
            } else if (col < col_to_remove) {
                reducedImage.at<Vec3b>(row, col) = image.at<Vec3b>(row, col);
            }
        }
    }

    return reducedImage;
}


/*
 * Remove seam from the gray image.
 */
Mat remove_seam_gray(Mat gimage, int seam[])
{
    int rows = gimage.rows;
    int cols = gimage.cols;
    Mat reducedImage(rows,cols-1,CV_8UC1);

    for (int row = 0; row < rows; row++) {
        int col_to_remove = seam[row];
        for (int col = 0; col < cols; col++) {
            if (col > col_to_remove) {
                reducedImage.at<uchar>(row, col-1) = gimage.at<uchar>(row, col);
            } else if (col < col_to_remove) {
                reducedImage.at<uchar>(row, col) = gimage.at<uchar>(row, col);
            }
        }
    }

    return reducedImage;
}


double min_helper(double x, double y, double z) {
    double min_val = x;
    if (y < min_val) {
        min_val = y;
    }
    if (z < min_val) {
        min_val = z;
    }
    return min_val;
}

/*
 * Find the seam of minimum energy according to our improved graph cut algorithm
 */
int *find_seam(Mat energy)
{
    cout << "finding seam" << endl;
    timestamp_t t0 = get_timestamp();

    int rows = energy.rows;
    int cols = energy.cols;
    int *seam = new int[rows];
    
    double min_val = -1;
    double min_col = -1;
    double inf = 100000;
    for (int col = 0; col < cols-1; col++) {
        double cur = energy.at<uchar>(rows-1, col);
        if (cur < min_val) {
            min_val = cur;
            min_col = col;
        }
    }
    int cur_col = min_col;
    int row;
    for (row = rows-1; row > 0; row--) {
        seam[row] = cur_col;
        double upleft, upright;
        if (cur_col != 0) {
            upleft = energy.at<uchar>(row-1, cur_col-1);
        } else {
            upleft = inf;
        }

        if (cur_col != cols-1) {
            upright = energy.at<uchar>(row-1, cur_col-1);
        } else {
            upright = inf;
        }

        double up = energy.at<uchar>(row-1, cur_col);

        min_val = min_helper(upleft, up, upright);

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
    return seam;
}

/*
 * Remove one vertical seam from img. The algorithm:
1) get energy matrix.
2) generate accumulated energy matrix
3) find seam
4) remove this seam from the image
 */
Mat reduce_vertical(Mat &grayImg1, Mat img) {
    int rows = grayImg1.rows;
    Mat energy = get_energy_matrix(grayImg1);
    Mat acm = generate_ACM(energy);
    int *seam = find_seam(acm);
    Mat returnImg = remove_seam(img, seam);
    grayImg1 = remove_seam_gray(grayImg1, seam);
    return returnImg;
}

/*
*  returns the energy matrix 
*/
Mat get_energy_matrix(Mat grayImg) {
    cout << "generating energy matrix" << endl;
    timestamp_t t0 = get_timestamp();
    int rows = grayImg.rows;
    int cols = grayImg.cols;
    Mat energy(rows, cols, CV_32S);
    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            double gradient = 0;
            if ((row == 0) | (row == rows - 1) | (cols == 0) | (col == cols -1)) {
                gradient = 1;
            } else {
                int right = grayImg.at<uchar>(row, col+1);
                int left = grayImg.at<uchar>(row, col-1);
                int up = grayImg.at<uchar>(row-1, col);
                int down = grayImg.at<uchar>(row+1, col);
                // right- left
                int dx = abs(right - left);
                // up - down
                int dy = abs(up - down);
                double diff = dx + dy;
                gradient = ((double)diff) / ((double)(255 * 2));
            }
            energy.at<uchar>(row, col) = gradient;
        }
    }
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "Energy matrix generation time: " << secs << "secs" << endl;
    return energy;
}

Mat generate_ACM(Mat energy) { 
    cout << "Generating ACM" << endl;
    timestamp_t t0 = get_timestamp();
    int rows = energy.rows;
    int cols = energy.cols;
    double inf = 100000;
    for (int row = 1; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            double up = energy.at<uchar>(row-1, col);
            if (col == 0) {
                double upright = energy.at<uchar>(row-1, col+1);
                energy.at<uchar>(row, col) = energy.at<uchar>(row, col) + min_helper(up, upright, inf);
            } else if (col == cols -1 ){
                double upleft = energy.at<uchar>(row-1, col-1);
                energy.at<uchar>(row, col) = energy.at<uchar>(row, col) + min_helper(up, upleft, inf);
            } else {
                double upright = energy.at<uchar>(row-1, col+1);
                double upleft = energy.at<uchar>(row-1, col-1);
                energy.at<uchar>(row, col) = energy.at<uchar>(row, col) + min_helper(up, upleft, upright);
            }
        }
    }
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "ACM generation time: " << secs << "secs" << endl;
    return energy;
}


/*
 * Remove one horizontal seam from img. The algorithm:
1) get energy matrix.
2) generate accumulated energy matrix
3) find seam
4) remove this seam from the image
 */
Mat reduce_horizontal(Mat &grayImg1, Mat img) {
    int rows = grayImg1.rows;
    Mat energy = get_energy_matrix(grayImg1.t());
    Mat acm = generate_ACM(energy);
    int *seam = find_seam(acm);
    Mat returnImg = remove_seam(img, seam);
    Mat grayImg1temp = remove_seam_gray(grayImg1.t(), seam);
    grayImg1 = grayImg1temp.t();
    return returnImg.t();
}

/*
 * Seam-cut to the image.
 * v: number of vertical seams to be removed
 * h: number of horizontal seams to be removed
 */
Mat reduce_frame(Mat image,  int v, int h) {
    Mat grayImg1;
    cvtColor(image, grayImg1, CV_RGB2GRAY);

    int min = 0, diff = 0;
    Mat reducedGrayImg1, reducedImg;
    
    grayImg1.copyTo(reducedGrayImg1);
    image.copyTo(reducedImg);
    
    if(h > v) {
        diff = h - v;
        min = v;
    } else {
        diff = v - h;
        min = h;
    }
    
    for(int i = 0; i < min; ++i) {
        reducedImg= reduce_vertical(reducedGrayImg1, reducedImg);
        reducedImg= reduce_horizontal(reducedGrayImg1,reducedImg.t());
    }
    
    if(h > v) {
        for(int i = 0; i < diff; ++i) {
            reducedImg= reduce_horizontal(reducedGrayImg1,reducedImg.t());
        }
    } else {
        for(int i = 0; i < diff; ++i) {
            reducedImg= reduce_vertical(reducedGrayImg1,reducedImg);
        }
    }
    return reducedImg;
}