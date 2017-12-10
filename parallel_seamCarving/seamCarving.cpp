#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cstring>
#include <algorithm>
#include <omp.h>
#include <limits.h>
#include "mic.h"

#define GET_INDEX(x,y,w) (x * w + y)*3
#define GRAY_SCALE(r,g,b) (0.21 * r + 0.71 * g + 0.08 * b) 
#define PI 3.14159265

using namespace std;

typedef struct {
  int index;
  int energy;
} seam_idx_t; 

static int _argc;
static const char **_argv;

using namespace std::chrono;
typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::duration<double> dsec;

extern void  png_to_text(const char *filename);
void calculate_energy(uint8_t *image, int *energy, int rows, int cols);
int reduce_image(uint8_t *image, int seam_count, int rows, int cols);
void remove_seam(uint8_t *outImage, uint8_t *image_temp, int *seams, seam_idx_t *seam_energy, int rows, int cols, int v, int seams_found);
int find_seams(int *energy, char *dir_map, int *seams, seam_idx_t *seam_energy, int rows, int cols);
bool bound_check(int row, int col, int rows, int cols);
void draw_seam(uint8_t *outImage, int *seams, seam_idx_t *seam_energy, int rows, int cols, int seams_found);
bool compare(const seam_idx_t &a, const seam_idx_t &b);


/* Starter code function, don't touch */
const char *get_option_string(const char *option_name,
                  const char *default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return _argv[i + 1];
  return default_value;
}

/* Starter code function, do not touch */
int get_option_int(const char *option_name, int default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return atoi(_argv[i + 1]);
  return default_value;
}

/* Starter code function, do not touch */
float get_option_float(const char *option_name, float default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return (float)atof(_argv[i + 1]);
  return default_value;
}

static void show_help(const char *program_path)
{
    printf("Usage: %s OPTIONS\n", program_path);
    printf("\n");
    printf("OPTIONS:\n");
    printf("\t-f <input_filename> (required)\n");
    printf("\t-w <width> (required)\n");
    printf("\t-h <height> (required)\n");
    printf("\t-s <num_of_seams> (required)\n");
    printf("\t-n <num_threads>\n");
}



void calculate_energy(uint8_t *image, char *dir_map, int *energy, int rows, int cols){

    int inf = 10000;
    #pragma omp parallel for
    for(int row = 0; row < rows; row++){
        for(int col = 0; col < cols; col++){
            int energy_val;
            int dir_val;
            
            if (row == rows - 1|| col == cols - 1 || col == 0 || row == 0){
                energy_val = inf;
                dir_val = 2;
            } else {

                // only looking at dx to improve performance
                // output same/better than when we looked at dx+dy


                //int rUp = image[GET_INDEX(row - 1, col, cols)];
                //int rDown = image[GET_INDEX(row + 1, col, cols)];
                int rLeft = image[GET_INDEX(row, col - 1, cols)];
                int rRight = image[GET_INDEX(row, col + 1, cols)];

                //int gUp = image[GET_INDEX(row - 1, col, cols) + 1];
                //int gDown = image[GET_INDEX(row + 1, col, cols) + 1];
                int gLeft = image[GET_INDEX(row, col - 1, cols) + 1];
                int gRight = image[GET_INDEX(row, col + 1, cols) + 1];

                //int bUp = image[GET_INDEX(row - 1, col, cols) + 2];
                //int bDown = image[GET_INDEX(row + 1, col, cols) + 2];
                int bLeft = image[GET_INDEX(row, col - 1, cols) + 2];
                int bRight = image[GET_INDEX(row, col + 1, cols) + 2];

                //double grayUp = GRAY_SCALE(rUp, gUp, bUp);
                //double grayDown = GRAY_SCALE(rDown, gDown, bDown);
                double grayRight = GRAY_SCALE(rRight, gRight, bRight);
                double grayLeft = GRAY_SCALE(rLeft, gLeft, bLeft);

                double gray_dx = grayRight - grayLeft;
                //double gray_dy = -(grayDown - grayUp);

                //double val = gray_dy/gray_dx;
                double val = gray_dx;
                double angle = static_cast<double>((atan(val) * 180/PI + 270)); 
                
                if (angle < 240) { 
                  dir_val = -1;
                } else if (angle > 300) {
                  dir_val = 1;
                } else {
                  dir_val = 0;
                }
                energy_val = abs(gray_dx);
                //energy_val = abs(gray_dx) + abs(gray_dy);
            }
            energy[row * cols + col] = energy_val;
            dir_map[row * cols +col] = dir_val;
        }
    }
}


bool bound_check(int row, int col, int rows, int cols) {
    if (col < 1 || col >= cols-1) return false;
    if (row < 1) return false;
    return true;
}

bool bound_check_2(int row, int col, int rows, int cols, int start, int end) {
    if (col < 1 || col >= cols-1) return false;
    if (col < start || col >= end) return false;
    if (row < 1) return false;
    return true;
}

/*
int find_seams(int *energy, char *dir_map, int *seams, seam_idx_t *seam_energy, int rows, int cols)
{   

    // need to reinitialize seams_energy list every time due to faulty random shuffle generator
    #pragma omp parallel for
    for (int i = 0; i < cols; i++) {
        seam_energy[i].index  = i;
        seam_energy[i].energy = INT_MAX;
    }

    //priority lists
    int plist_straight[3] {0, 1, -1};
    int plist_left[3] {-1, 0, 1};
    int plist_right[3] {1, 0, -1};
    
    int num_priorities = 3;

    int rand_order[cols-1];
    #pragma omp parallel for
    for (int i = 1; i < cols; i++){
        rand_order[i] = i;
    }

    // randomize the start cols
    random_shuffle(&rand_order[0], &rand_order[cols-1]);

    // the random_shuffle seams to give broken output at times
    // this check is there to make sure a seg fault does not occur
    #pragma omp parallel for
    for (int i = 1; i < cols; i++){
        if (rand_order[i] > cols-1) {
            rand_order[i] = -1;
            printf("Incorrect random number found.\n");
        }
    }

    int count = 0;

    int seam_col, col, cur_energy, row, dir;

    int NITEMS = rows*cols;
    omp_lock_t lock[NITEMS];

    #pragma omp parallel for
    for (int i=0; i<NITEMS; i++)
        omp_init_lock(&(lock[i]));

    #pragma omp parallel for private(seam_col, col, row, dir) 
    for (int i = 1; i < cols; i+=1){

        seam_col = rand_order[i];
        if (seam_col > cols-1) {
          //ignore this one
          continue;
        }
        
        col = seam_col; 
        
        dir = dir_map[col];
        seams[seam_col * rows] = col;

        cur_energy = 0;

        for (row = 0; row < rows - 2; row++){

            cur_energy += energy[row * cols + col];

            int temp_col;
            
            if (dir == 0) {
                bool found = false;
                //priority list of len 3
                for (int i = 0; i < num_priorities; i++) {
                  temp_col = col + plist_straight[i];
                  if (bound_check((row+1), temp_col, rows, cols)) {
                    int index = (row+1) * cols + temp_col;
                    
                    //check if already used by another seam
                    if (dir_map[index] != 2) {
                        omp_set_lock(&(lock[index]));
                        if (dir_map[index] != 2) {
                            dir = dir_map[index];
                            dir_map[index] = 2;
                            omp_unset_lock(&(lock[index]));
                            found = true;
                            col = temp_col;
                            seams[seam_col * rows + (row+1)] = col;
                            break;
                        }
                        omp_unset_lock(&(lock[index]));
                    }    
                  }
                }
                if (!found) {
                  cur_energy = -1;
                  break;
                }

            }else if(dir == -1) {

                bool found = false;
                //priority list of len 3
                for (int i = 0; i < num_priorities; i++) {
                  temp_col = col + plist_left[i];
                  if (bound_check((row+1), temp_col, rows, cols)) {
                    int index = (row+1) * cols + temp_col;
                    
                    //check if already used by another seam
                    if (dir_map[index] != 2) {
                        omp_set_lock(&(lock[index]));
                        if (dir_map[index] != 2) {
                            dir = dir_map[index];
                            dir_map[index] = 2;
                            omp_unset_lock(&(lock[index]));
                            found = true;
                            col = temp_col;
                            seams[seam_col * rows + (row+1)] = col;
                            break;
                        }
                        omp_unset_lock(&(lock[index]));
                    }    
                  }
                }
                if (!found) {
                  cur_energy = -1;
                  break;
                }

            } else {

              // dir = 1
                bool found = false;
                //priority list of len 3
                for (int i = 0; i < num_priorities; i++) {
                  temp_col = col + plist_right[i];
                  if (bound_check((row+1), temp_col, rows, cols)) {
                    int index = (row+1) * cols + temp_col;
                    
                    //check if already used by another seam
                    if (dir_map[index] != 2) {
                        omp_set_lock(&(lock[index]));
                        if (dir_map[index] != 2) {
                            dir = dir_map[index];
                            dir_map[index] = 2;
                            omp_unset_lock(&(lock[index]));
                            found = true;
                            col = temp_col;
                            seams[seam_col * rows + (row+1)] = col;
                            break;
                        }
                        omp_unset_lock(&(lock[index]));
                    }    
                  }
                }
                if (!found) {
                  cur_energy = -1;
                  break;
                }
            }
        }

        if (cur_energy != -1) {
            seam_energy[seam_col].index = seam_col;
            seam_energy[seam_col].energy = cur_energy;
            #pragma omp atomic 
            count++;
        } 
    }

    #pragma omp parallel for
    for (int i=0; i<NITEMS; i++)
        omp_destroy_lock(&(lock[i]));

    return count;
}
*/


int find_seams(int *energy, char *dir_map, int *seams, seam_idx_t *seam_energy, int rows, int cols)
{   

    //priority lists
    int plist_straight[3] {0, 1, -1};
    int plist_left[3] {-1, 0, 1};
    int plist_right[3] {1, 0, -1};
    
    int num_priorities = 3;

    int count = 0;

    int num_threads = 16;
    int offset = cols / num_threads;
 
    #pragma omp parallel for  
    for (int thread = 0; thread < num_threads; thread++) 
    {
        int seam_col, col, cur_energy, row, dir;
        int start = offset * thread;
        int end;
        if (thread !=  num_threads -1){
          end = offset *  (thread+1);
        } else {
          end = cols;
        }
      
        for (int i = start; i < end; i+=1) {

          seam_col = i;
          col = seam_col;  
          dir = dir_map[col];
          seams[seam_col * rows] = col;
          cur_energy = 0;

            for (row = 0; row < rows - 2; row++){

                cur_energy += energy[row * cols + col];

                int temp_col;
                
                if (dir == 0) {
                    bool found = false;
                    //priority list of len 3
                    for (int i = 0; i < num_priorities; i++) {
                      temp_col = col + plist_straight[i];
                      if (bound_check_2((row+1), temp_col, rows, cols, start, end)) {
                        int index = (row+1) * cols + temp_col;
                        
                        //check if already used by another seam
                        if (dir_map[index] != 2) {
                            dir = dir_map[index];
                            dir_map[index] = 2;
                            found = true;
                            col = temp_col;
                            seams[seam_col * rows + (row+1)] = col;
                            break;
                        }    
                      }
                    }
                    if (!found) {
                      cur_energy = -1;
                      break;
                    }

                }else if(dir == -1) {

                    bool found = false;
                    //priority list of len 3
                    for (int i = 0; i < num_priorities; i++) {
                      temp_col = col + plist_left[i];
                      if (bound_check_2((row+1), temp_col, rows, cols, start, end)) {
                        int index = (row+1) * cols + temp_col;
                        
                        //check if already used by another seam
                        if (dir_map[index] != 2) {
                            dir = dir_map[index];
                            dir_map[index] = 2;
                            found = true;
                            col = temp_col;
                            seams[seam_col * rows + (row+1)] = col;
                            break;
                        }   
                      }
                    }
                    if (!found) {
                      cur_energy = -1;
                      break;
                    }

                } else {

                  // dir = 1
                    bool found = false;
                    //priority list of len 3
                    for (int i = 0; i < num_priorities; i++) {
                      temp_col = col + plist_right[i];
                      if (bound_check_2((row+1), temp_col, rows, cols, start, end)) {
                        int index = (row+1) * cols + temp_col;
                        
                        //check if already used by another seam
                        if (dir_map[index] != 2) {
                            dir = dir_map[index];
                            dir_map[index] = 2;
                            found = true;
                            col = temp_col;
                            seams[seam_col * rows + (row+1)] = col;
                            break;
                        }    
                      }
                    }
                    if (!found) {
                      cur_energy = -1;
                      break;
                    }
                }
            }

          if (cur_energy != -1) {
            seam_energy[seam_col].index = seam_col;
            seam_energy[seam_col].energy = cur_energy;
            #pragma omp atomic 
            count++;
          } else {
            seam_energy[seam_col].index  = seam_col;
            seam_energy[seam_col].energy = INT_MAX;
          }
        }
    }
    return count;
}



bool compare(const seam_idx_t &a, const seam_idx_t &b)
{
    return a.energy < b.energy;
}

/*
 * Remove seam from the image.
 */
void remove_seam(uint8_t *image, uint8_t *image_temp, int *seams, seam_idx_t *seam_energy, int rows, int cols, int batch_size, int seams_found)
{   


    sort(seam_energy, seam_energy + (cols), compare);

    char temp[rows*cols] = {0};

    int j=0;
    for (int i = 0; i < seams_found && i < batch_size; i++) {
        j++;
        int start_col = seam_energy[i].index;

        #pragma omp parallel for
        for (int row = 0; row < rows; row++) {

            int col_to_remove = seams[start_col * rows + row];

            for (int col = 0; col < cols; col++) {

                int index = row * cols + col;
                if (col == col_to_remove) {

                    temp[index] = 1;
                    j++;

                } else {

                    if (temp[index] != 1) temp[index] = 0;

                }
            }
        }
    } 
    int count = 0;
    #pragma omp parallel for
    for (int row = 0; row < rows; row++){

         int offset = 0;

        for (int col = 0; col < cols; col++) {

            int index = row * cols + col;

            if (temp[index] == 1) {

                offset++;
                count++;
            } else {

                int index = row * cols + col;
                int index3 = index*3;
                int new_index = row * (cols - batch_size) + col-offset;
                int new_index3 = new_index * 3; 

                image_temp[new_index3] = image[index3];
                image_temp[new_index3+1] = image[index3+1];
                image_temp[new_index3+2] = image[index3+2];
            }
        }
    }

    int new_width = cols-batch_size;

    memcpy(image, image_temp, sizeof(uint8_t) * new_width * rows * 3);
}



/*
 * Draw seam on the image.
 */
void draw_seam(uint8_t *outImage, int *seams, seam_idx_t *seam_energy, int rows, int cols, int seams_found)
{

    sort(seam_energy, seam_energy + (cols), compare);
     
     for (int i = 0; i < cols; i++) {

        if (seam_energy[i].energy == INT_MAX) continue;

        int start_col = seam_energy[i].index;

        for (int row = 0; row < rows; row++) {

            int col_to_remove = seams[start_col * rows + row];

            for (int col = 0; col < cols; col++) {
                int index = row * cols + col;
                int index3 = index*3;
                if (col == col_to_remove) {
                    outImage[index3] = 255;
                    outImage[index3 + 1] = 0;
                    outImage[index3 + 2] = 0;
                }
            }
        }
    } 
}


int reduce_image(uint8_t *reducedImg, uint8_t *image_temp, int *energy, int *seam, int v, int rows, int cols) {

    
   double energy_compute_time = 0;
   double seam_finding_time = 0;
   double seam_removal_time = 0;

   int batch_size = 30;
   char *dir_map = (char *)calloc(rows * cols, sizeof(char));
   int *seams = (int *)calloc(rows * cols, sizeof(int));
   seam_idx_t *seam_energy = (seam_idx_t *)calloc(cols, sizeof(seam_idx_t));


   for (int i = 0; i < v; i += batch_size) {

      int cur_size; 

      if (i + batch_size > v) {
          cur_size = v - i;
      } else {
          cur_size = batch_size;
      }
      

      auto now = Clock::now();    
      calculate_energy(reducedImg, dir_map, energy, rows, cols);
      energy_compute_time += duration_cast<dsec>(Clock::now() - now).count();

      now = Clock::now();
      int seams_found = find_seams(energy, dir_map, seams, seam_energy, rows, cols);
      seam_finding_time += duration_cast<dsec>(Clock::now() - now).count();

      printf("seams found: %d\n", seams_found);
      if (seams_found < cur_size) {
        printf("Error! Not enough seams found\n");
        return cols;
      }

      now = Clock::now();
      remove_seam(reducedImg, image_temp, seams, seam_energy, rows, cols, cur_size, seams_found);
      seam_removal_time += duration_cast<dsec>(Clock::now() - now).count();

      cols -= cur_size;
   }
/*
    auto now = Clock::now();    
    calculate_energy(reducedImg, dir_map, energy, rows, cols);
    energy_compute_time += duration_cast<dsec>(Clock::now() - now).count();

    now = Clock::now();
    int seams_found = find_seams(energy, dir_map, seams, seam_energy, rows, cols);
    seam_finding_time += duration_cast<dsec>(Clock::now() - now).count();

    printf("seams found: %d\n", seams_found);
    now = Clock::now();
    draw_seam(reducedImg, seams, seam_energy, rows, cols, seams_found); 
    seam_removal_time += duration_cast<dsec>(Clock::now() - now).count();
*/ 
   printf("Total energy calculation Time: %lf.\n", energy_compute_time);
   printf("Total seam finding Time: %lf.\n", seam_finding_time);
   printf("Total seam removal Time: %lf.\n", seam_removal_time);
   return cols;
}


int main(int argc, const char *argv[])
{
  auto init_start = Clock::now();
  double init_time = 0;
 
  _argc = argc - 1;
  _argv = argv + 1;

  /* You'll want to use these parameters in your algorithm */
  const char *input_filename = get_option_string("-f", NULL);
  int seam_count = get_option_int("-s", -1);
  int num_of_threads = get_option_int("-n", 1);
  int error = 0;

  if (input_filename == NULL) {
    printf("Error: You need to specify -f.\n");
    error = 1;
  }

  if (seam_count == -1) {
    printf("Error: You need to specify -s.\n");
    error = 1;
  }

  if (error) {
    show_help(argv[0]);
    return 1;
  }

  // Set our thread count.
  omp_set_num_threads(num_of_threads);

  printf("Number of threads: %d\n", num_of_threads);

  FILE *input = fopen(input_filename, "r");

  if (!input) {
    printf("Unable to open file: %s.\n", input_filename);
    return -1;
  }

  int r, g, b;
  int index = 0;
  int rows, cols;
  fscanf(input, "%d %d\n", &cols, &rows);

  int img_size = rows * cols * 3;

  uint8_t *image = (uint8_t *)calloc(img_size, sizeof(uint8_t));
  uint8_t *image_temp = (uint8_t *)calloc(img_size, sizeof(uint8_t));
  int *energy = (int *)calloc(cols * rows, sizeof(int));
  int *seam = (int *)calloc(rows, sizeof(int));

  while (fscanf(input, "%d %d %d\n", &r, &g, &b) != EOF) {
    /* PARSE THE INPUT FILE HERE */
    image[index] = (uint8_t)r; index++;
    image[index] = (uint8_t)g; index++;
    image[index] = (uint8_t)b; index++;
  }
  fclose(input);
  init_time += duration_cast<dsec>(Clock::now() - init_start).count();
  printf("Initialization Time: %lf.\n", init_time);
  auto compute_start = Clock::now();
  double compute_time = 0;

#ifdef RUN_MIC /* Use RUN_MIC to distinguish between the target of compilation */

  /* This pragma means we want the code in the following block be executed in 
   * Xeon Phi.
   */
#pragma offload target(mic) \
  inout(image: length(img_size) INOUT) \
  inout(image_temp: length(img_size) INOUT) \
  inout(energy: length(rows * cols) INOUT) \
  inout(seam: length(rows) INOUT) 
#endif
  {
    /* Implement the wire routing algorithm here
     * Feel free to structure the algorithm into different functions
     * Don't use global variables.
     * Use OpenMP to parallelize the algorithm. 
     * You should really implement as much of this (if not all of it) in
     * helper functions. */
     //initialize_costs(costs, wires, num_of_wires, dim_x, dim_y);
    //int seams_removed = reduce_image(image, image_temp, energy, seam, seam_count, rows, cols);
  }
  int new_width = reduce_image(image, image_temp, energy, seam, seam_count, rows, cols);
  compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
  printf("Computation Time: %lf.\n", compute_time);
  printf("Total time: %1f.\n", compute_time + init_time);
  
  /* OUTPUT YOUR RESULTS TO FILES HERE 
   * When you're ready to output your data to files, uncommment this chunk of
   * code and fill in the specified blanks indicated by comments. More about
   * this in the README. */  

  FILE *outFile = fopen("outputImg.txt", "w");
  if (!outFile) {
    printf("Error: couldn't output image file");
    return -1;
  }
  // output information here
  fprintf(outFile,"%d %d\n", new_width, rows);
  
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < new_width; col++) {
      //printf("%d/n");
      int index = row * new_width + col;
      int index3 = index * 3;
      fprintf(outFile,"%d %d %d\n", (int)image[index3], (int)image[index3+1], (int)image[index3+2]);
    }
  }
  free(image);
  fclose(outFile);

  return 0;
}