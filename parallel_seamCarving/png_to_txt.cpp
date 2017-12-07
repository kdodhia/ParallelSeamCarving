#include <stdlib.h>
#include <stdio.h>
#include <png.h>
#include <stdint.h>


int main(int argc, const char *argv[])
{
  const char *filename = "image.png"
  FILE *outputImageFile = fopen("image_rgb.txt", "w");
  if (!outputImageFile) return -1;

  FILE *fp = fopen(filename, "rb");

  png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  
  if(!png) abort();

  png_infop info = png_create_info_struct(png);
  if(!info) abort();

  if(setjmp(png_jmpbuf(png))) abort();

  png_init_io(png, fp);

  png_read_info(png, info);

  int width = png_get_image_width(png, info);
  int height = png_get_image_height(png, info);
  png_byte color_type = png_get_color_type(png, info);
  png_byte bit_depth  = png_get_bit_depth(png, info);
  (void)color_type;
  (void)bit_depth;

  printf("width:%d\n", width);
  printf("height:%d\n", height);

  uint8_t *r = (uint8_t*)calloc(width * height, sizeof(uint8_t));
  uint8_t *g = (uint8_t*)calloc(width * height, sizeof(uint8_t));
  uint8_t *b = (uint8_t*)calloc(width * height, sizeof(uint8_t));

  png_bytep row_pointers[height];

  for(int y = 0; y < height; y++) {
    row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
  }

  png_read_image(png, row_pointers);

  for (int row = 0; row < height; row++) {

     // Get the pointer for this row
     png_byte *thisRowPointer = (png_byte*)row_pointers[row];

     for (int col = 0; col < width; col++) {
          int colIndex = 3 * col;
          r[(row * width) + col] = thisRowPointer[colIndex];
          g[(row * width) + col] = thisRowPointer[colIndex + 1];
          b[(row * width) + col] = thisRowPointer[colIndex + 2];
     }
  }

  fprintf(outputImageFile,"%d %d\n", width, height);
  for (int i = 0; i < width * height; i++) {
   fprintf(outputImageFile,"%d %d %d\n", r[i], g[i], b[i]);
  }

  for(int y = 0; y < height; y++) {
    free(row_pointers[y]);
  }

  fclose(fp);
  fclose(outputImageFile);
}
