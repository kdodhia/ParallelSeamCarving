#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>

#include <pthread.h>

#include "graph.h"

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

    string::size_type pAt = inFile.find_last_of('.');   // Find extension point
    const string outFile = inFile.substr(0, pAt) + "-result.jpg";
    outImage = reduce();
    imwrite(outFile, outImage);
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "Total runtime: " << secs << "secs" << endl;
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
 *
 * image: image to be seam-removed
 * seam: column index of each row
 */
Mat remove_seam(Mat image, int seam[])
{
    int nrows = image.rows;
    int ncols = image.cols;
    Mat reducedImage(nrows,ncols-1,CV_8UC3);
    
    for(int i=0; i<nrows; i++)
    {
        if(seam[i] != 0)
            image.row(i).colRange(Range(0,seam[i])).copyTo(reducedImage.row(i).colRange(Range(0,seam[i])));
        if(seam[i] != ncols-1)
            image.row(i).colRange(Range(seam[i]+1, ncols)).copyTo(reducedImage.row(i).colRange(Range(seam[i],ncols-1)));
    }
    return reducedImage;
}


/*
 * Remove seam from the gray image.
 *
 * GrayImage: gray image to be seam-removed
 * seam: column index of each row
 */
Mat remove_seam_gray(Mat GrayImage, int seam[])
{
    int nrows = GrayImage.rows;
    int ncols = GrayImage.cols;
    Mat reducedImage(nrows,ncols-1,CV_8UC1);
    for(int i=0; i<nrows; i++)
    {
        if(seam[i] != 0)
            GrayImage.row(i).colRange(Range(0,seam[i])).copyTo(reducedImage.row(i).colRange(Range(0,seam[i])));
        if(seam[i] != ncols-1)
            GrayImage.row(i).colRange(Range(seam[i]+1, ncols)).copyTo(reducedImage.row(i).colRange(Range(seam[i],ncols-1)));
    }
    return reducedImage;
}


/*
 * Find the seam of minimum energy according to our improved graph cut algorithm
 *
 * grayImg1: current image needed to be seam-removed
 */
int *find_seam(Mat grayImg1)
{
    cout << "finding seam" << endl;
    timestamp_t t0 = get_timestamp();
    // define variables
    typedef Graph<int,int,int> GraphType;
    int rows = grayImg1.rows;
    int cols = grayImg1.cols;
    double inf = 100000;
    int *Seam = new int[rows];
    float a1= 1;
    GraphType *g = new GraphType(/*estimated # of nodes*/ rows*cols, /*estimated # of edges*/ ((rows-1)*cols + (cols-1)*rows + 2*(rows-1)*(cols-1)));
    
    int LR, LR1;
    int posLU, posLU1;
    int negLU, negLU1;
    for (int i = 1; i<=rows*cols; i++) {
        g -> add_node();
    }
    
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            if(j==0) {
                g -> add_tweights( i*cols,   /* capacities */  inf,0 );
            }
            else if(j==cols-1) {
                g -> add_tweights( ((i+1)*cols) -1,   /* capacities */ 0, inf );
            }
            
            if(j==0) {
                LR1= grayImg1.at<unsigned char>(i,j+1);
                LR= a1*LR1;
                g -> add_edge( i*cols, i*cols+1,    /* capacities */ LR, inf );
            }
            else if(j!=cols-1) {
                LR1= abs(grayImg1.at<unsigned char>(i,j+1)-grayImg1.at<unsigned char>(i,j-1));
                LR= a1*LR1;
                g -> add_edge( i*cols + j, i*cols + j +1, LR, inf );
            }
            
            if(i!=rows-1) {
                if(j==0) {
                    // positive LU
                    posLU1= grayImg1.at<unsigned char>(i,j);
                    posLU= a1*posLU1;
                    // negative LU
                    negLU1= grayImg1.at<unsigned char>(i+1,j);
                    negLU= a1*negLU1;
                    g -> add_edge( i*cols + j, i*cols + j +1, negLU, posLU );
                }
                else {
                    // positive LU
                    posLU1= abs(grayImg1.at<unsigned char>(i,j)-grayImg1.at<unsigned char>(i+1,j-1));
                    posLU= a1*posLU1;
                    // negative LU
                    negLU1= abs(grayImg1.at<unsigned char>(i+1,j)-grayImg1.at<unsigned char>(i,j-1));
                    negLU= a1*negLU1;
                    g -> add_edge( i*cols + j, i*cols + j +1, negLU, posLU );
                }
            }
            if(i!=0 && j!=0)
            {
                g -> add_edge( i*cols + j, (i-1)*cols + j-1, inf, 0 );
            }
            if(i!=rows-1 && j!=0)
            {
                g -> add_edge( i*cols + j, (i+1)*cols + j-1, inf, 0 );
            }
        }
    }
    
    int flow = g -> maxflow();
    for(int i=0; i<rows; i++) {
        for(int j=0;j<cols; j++) {
            if(g->what_segment(i*cols+j) == GraphType::SINK) {
                Seam[i] = j-1;
                break;
            }
            if(j==cols-1 && g->what_segment(i*cols+j) == GraphType::SOURCE) {
                Seam[i] = cols-1;
            }
        }
    }
    delete g;
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    cout << "Seam generation time: " << secs << "secs" << endl;
    return Seam;
}


/*
 * Remove one vertical seam from img.
 * 
 * grayImg1: used to calculate the seam
 * img: image needed to perform seam-cut
 */
Mat reduce_vertical(Mat &grayImg1, Mat img) {
    int rows = grayImg1.rows;
    int *seam = new int[rows];
    seam = find_seam(grayImg1);
    Mat returnImg = remove_seam(img, seam);
    grayImg1 = remove_seam_gray(grayImg1, seam);
    return returnImg;
}


/*
 * Remove one horizontal seam from img.
 *
 * grayImg1: used to calculate the seam
 * img: image needed to perform seam-cut
 */
Mat reduce_horizontal(Mat &grayImg1, Mat img) {
    int rows = grayImg1.rows;
    int *seam = new int[rows];
    seam = find_seam(grayImg1.t());
    Mat returnImg = remove_seam(img, seam);
    Mat grayImg1temp = remove_seam_gray(grayImg1.t(), seam);
    grayImg1 = grayImg1temp.t();
    return returnImg.t();
}

/*
 * Seam-cut to the frame1.
 *
 * frame1: the frame to be seam-cut
 * v: number of vertical seams to be removed
 * h: number of horizontal seams to be removed
 */
Mat reduce_frame(Mat frame1,  int v, int h) {
    Mat grayImg1;
    cvtColor(frame1,grayImg1, CV_RGB2GRAY);

    int min = 0, diff = 0;
    Mat reducedGrayImg1, reducedImg;
    
    grayImg1.copyTo(reducedGrayImg1);
    frame1.copyTo(reducedImg);
    
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

