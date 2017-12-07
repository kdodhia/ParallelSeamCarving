#include <stdlib.h>
#include <stdio.h>
#include <png.h>
#include <stdint.h>


int main(int argc, const char *argv[])
{
    const char *imageName = "outputImg.txt";
    FILE *inputImageFile = fopen(imageName, "r");
    int width, height;
    fscanf(inputImageFile, "%d %d\n", width, height);

    printf("Width: %d Height: %d\n", width, height);

    uint8_t *r = (uint8_t *)calloc(width * height, sizeof(uint8_t));
    uint8_t *g = (uint8_t *)calloc(width * height, sizeof(uint8_t));
    uint8_t *b = (uint8_t *)calloc(width * height, sizeof(uint8_t));
    int rval, gval, bval;
    int i = 0;
    while (fscanf(inputImageFile, "%d %d %d\n", &rval, &gval, &bval) != EOF) {
      r[i] = (uint8_t)rval;
      g[i] = (uint8_t)gval;
      b[i] = (uint8_t)bval;
      i++;
    }

     png_structp png = NULL;
     png_infop info = NULL;
     png_byte *rowBytes;//[image.width * 3];
     FILE *fp = NULL;

     fp = fopen("output.png", "wb");
     if (fp == NULL) {
          printf("Error writing PNG image\n");
          return;
     }


     png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
     if (png == NULL) {
          printf("Error writing PNG image\n");
          return;
     }

     info = png_create_info_struct(png);
     if (info == NULL) {
          printf("Error writing PNG image\n");
          return;
     }

     png_init_io(png, fp);

     png_set_IHDR(png, info, width, height,
               8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

     png_text title_text;
     title_text.compression = PNG_TEXT_COMPRESSION_NONE;
     char key[] = "Title";
     char text[] = "output.png";
     title_text.key = key;
     title_text.text = text;
     png_set_text(png, info, &title_text, 1);

     png_write_info(png, info);


     rowBytes = (png_byte*)malloc(3 * width * sizeof(png_byte));

     for (int row = 0; row < image.height; row++) {
          for (int col = 0; col < image.width; col++) {
               // int colIndex = col * 3;
               rowBytes[(col * 3)]   = (png_byte)r[(width * row) + col];
               rowBytes[(col * 3) + 1] = (png_byte)g[(width * row) + col];
               rowBytes[(col * 3) + 2] = (png_byte)b[(width * row) + col];
          }
          png_write_row(png, (png_bytep)rowBytes);
     }

     png_write_end(png, NULL);
}