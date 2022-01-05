#include <stdbool.h>
#include "writeBMP.h"
#include "readBMP.h"
#include <stdlib.h>
#define KERNEL_SIZE 3


typedef struct {
    short red;
    short green;
    short blue;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;


/* Compute min and max of two integers, respectively */
int min(int a, int b) { return (a < b ? a : b); }
int max(int a, int b) { return (a > b ? a : b); }

int calcIndex(int i, int j, int n) {
	return ((i)*(n)+(j));
}

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
void initialize_pixel_sum(pixel_sum *sum) {
	sum->red = sum->green = sum->blue = 0;
	// sum->num = 0;
	return;
}

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

	// divide by kernel's weight
	sum.red = sum.red / kernelScale;
	sum.green = sum.green / kernelScale;
	sum.blue = sum.blue / kernelScale;

	// truncate each pixel's color values to match the range [0,255]
	current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
	current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
	current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
	return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {
	sum->red += ((int) p.red) * weight;
	sum->green += ((int) p.green) * weight;
	sum->blue += ((int) p.blue) * weight;
	// sum->num++;
	return;
}

/*
 *  Applies kernel for pixel at (i,j)
 */
static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	int ii, jj;
	int currRow, currCol;
	pixel_sum sum;
	pixel current_pixel;
	int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
	int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
	int min_row, min_col, max_row, max_col;
	pixel loop_pixel;

	initialize_pixel_sum(&sum);

	for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) {
		for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) {

			int kRow, kCol;

			// compute row index in kernel
			if (ii < i) {
				kRow = 0;
			} else if (ii > i) {
				kRow = 2;
			} else {
				kRow = 1;
			}

			// compute column index in kernel
			if (jj < j) {
				kCol = 0;
			} else if (jj > j) {
				kCol = 2;
			} else {
				kCol = 1;
			}

			// apply kernel on pixel at [ii,jj]
			sum_pixels_by_weight(&sum, src[calcIndex(ii, jj, dim)], kernel[kRow][kCol]);
		}
	}

	if (filter) {
		// find min and max coordinates
		for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) {
			for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) {
				// check if smaller than min or higher than max and update
				loop_pixel = src[calcIndex(ii, jj, dim)];
				if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) <= min_intensity) {
					min_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
					min_row = ii;
					min_col = jj;
				}
				if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) > max_intensity) {
					max_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
					max_row = ii;
					max_col = jj;
				}
			}
		}
		// filter out min and max
		sum_pixels_by_weight(&sum, src[calcIndex(min_row, min_col, dim)], -1);
		sum_pixels_by_weight(&sum, src[calcIndex(max_row, max_col, dim)], -1);
	}

	// assign kernel's result to pixel at [i,j]
	assign_sum_to_pixel(&current_pixel, sum, kernelScale);
	return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	int i, j;
	for (i = 1 ; i < dim - 1; i++) {
		for (j =  1 ; j < dim - 1 ; j++) {
			dst[calcIndex(i, j, dim)] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
		}
	}
}

void blurKernelSmoothWithoutFilter9(int dim, pixel* src, pixel* dst){
    int row, col;
    int b, r, g;

    //row 0
    for(col = 1; col < dim - 1; ++col){
        //initialize also
        dst[dim + col].blue = (src[col - 1].blue + src[col].blue + src[col + 1].blue);
        dst[dim + col].red = (src[col - 1].red + src[col].red + src[col + 1].red);
        dst[dim + col].green = (src[col - 1].green + src[col].green + src[col + 1].green);
    }

    //row 1
    for(col = 1; col < dim - 1; ++col){
        b = (src[dim + col -1].blue + src[dim + col].blue + src[dim + col + 1].blue);
        r = (src[dim + col -1].red + src[dim + col].red + src[dim + col + 1].red);
        g = (src[dim + col -1].green + src[dim + col].green + src[dim + col + 1].green);
        //add the value
        dst[dim + col].blue += b;
        dst[dim + col].red += r;
        dst[dim + col].green += g;
        //initialize also
        dst[2*dim + col].blue = b;
        dst[2*dim + col].red = r;
        dst[2*dim + col].green = g;
    }

    //rows 2 - (dim-3)
    for(row = 2; row <= dim - 3; ++row){
        for(col = 1; col < dim - 1; ++col){
          b = (src[row*dim + col - 1].blue + src[row*dim + col].blue + src[row*dim + col + 1].blue);
          g = (src[row*dim + col - 1].green + src[row*dim + col].green + src[row*dim + col + 1].green);
          r = (src[row*dim + col - 1].red + src[row*dim + col].red + src[row*dim + col + 1].red);
          //add the value
          dst[(row-1)*dim + col].blue += b;
          dst[(row-1)*dim + col].red += r;
          dst[(row-1)*dim + col].green += g;
          //close
          dst[(row-1)*dim + col].blue /= 9;
          dst[(row-1)*dim + col].red /= 9;
          dst[(row-1)*dim + col].green /= 9;
          //add the value
          dst[row*dim + col].blue += b;
          dst[row*dim + col].red += r;
          dst[row*dim + col].green += g;
          //initialize also
          dst[(row+1)*dim + col].blue = b;
          dst[(row+1)*dim + col].red = r;
          dst[(row+1)*dim + col].green = g;
        }
    }
    //row = dim -2;
    for(col = 1; col < dim -1; ++col){
        b = (src[(dim - 2)*dim + col -1].blue + src[(dim - 2)*dim + col].blue + src[(dim - 2)*dim + col + 1].blue);
        r = (src[(dim - 2)*dim + col -1].red + src[(dim - 2)*dim + col].red + src[(dim - 2)*dim + col + 1].red);
        g = (src[(dim - 2)*dim + col -1].green + src[(dim - 2)*dim + col].green + src[(dim - 2)*dim + col + 1].green);
        //add the value
        dst[(dim - 3)*dim + col].blue += b;
        dst[(dim - 3)*dim + col].red += r;
        dst[(dim - 3)*dim + col].green += g;
        //close
        dst[(dim - 3)*dim + col].blue /= 9;
        dst[(dim - 3)*dim + col].red /= 9;
        dst[(dim - 3)*dim + col].green /= 9;
        //add the value
        dst[(dim - 2)*dim + col].blue += b;
        dst[(dim - 2)*dim + col].red += r;
        dst[(dim - 2)*dim + col].green += g;
    }
    //row = dim - 1
    for(col = 1; col < dim -1; ++col){
        dst[(dim - 2)*dim + col].blue += (src[(dim - 1)*dim + col - 1].blue + src[(dim - 1)*dim + col].blue + src[(dim - 1)*dim + col + 1].blue);
        dst[(dim - 2)*dim + col].red += (src[(dim - 1)*dim + col - 1].red + src[(dim - 1)*dim + col].red + src[(dim - 1)*dim + col + 1].red);
        dst[(dim - 2)*dim + col].green += (src[(dim - 1)*dim + col - 1].green + src[(dim - 1)*dim + col].green + src[(dim - 1)*dim + col + 1].green);
        //close
        dst[(dim - 2)*dim + col].blue /= 9;
        dst[(dim - 2)*dim + col].red /= 9;
        dst[(dim - 2)*dim + col].green /= 9;
    }
}

void blurKernelSmoothWithFilter7(int dim, pixel* src, pixel* dst) {
    int row, col;
    int b, r, g;

    //row 0
    for (col = 1; col < dim - 1; ++col) {
        //initialize also
        dst[dim + col].blue = (src[col - 1].blue + src[col].blue + src[col + 1].blue);
        dst[dim + col].red = (src[col - 1].red + src[col].red + src[col + 1].red);
        dst[dim + col].green = (src[col - 1].green + src[col].green + src[col + 1].green);
    }

    //row 1
    for (col = 1; col < dim - 1; ++col) {
        b = (src[dim + col - 1].blue + src[dim + col].blue + src[dim + col + 1].blue);
        r = (src[dim + col - 1].red + src[dim + col].red + src[dim + col + 1].red);
        g = (src[dim + col - 1].green + src[dim + col].green + src[dim + col + 1].green);
        //add the value
        dst[dim + col].blue += b;
        dst[dim + col].red += r;
        dst[dim + col].green += g;
        //initialize also
        dst[2 * dim + col].blue = b;
        dst[2 * dim + col].red = r;
        dst[2 * dim + col].green = g;
    }

    //rows 2 - (dim-3)
    for (row = 2; row <= dim - 3; ++row) {
        for (col = 1; col < dim - 1; ++col) {
            b = (src[row * dim + col - 1].blue + src[row * dim + col].blue + src[row * dim + col + 1].blue);
            g = (src[row * dim + col - 1].green + src[row * dim + col].green + src[row * dim + col + 1].green);
            r = (src[row * dim + col - 1].red + src[row * dim + col].red + src[row * dim + col + 1].red);
            //add the value
            dst[(row - 1) * dim + col].blue += b;
            dst[(row - 1) * dim + col].red += r;
            dst[(row - 1) * dim + col].green += g;
            //add the value
            dst[row * dim + col].blue += b;
            dst[row * dim + col].red += r;
            dst[row * dim + col].green += g;
            //initialize also
            dst[(row + 1) * dim + col].blue = b;
            dst[(row + 1) * dim + col].red = r;
            dst[(row + 1) * dim + col].green = g;
        }
    }
    //row = dim -2;
    for (col = 1; col < dim - 1; ++col) {
        b = (src[(dim - 2) * dim + col - 1].blue + src[(dim - 2) * dim + col].blue + src[(dim - 2) * dim + col + 1].blue);
        r = (src[(dim - 2) * dim + col - 1].red + src[(dim - 2) * dim + col].red + src[(dim - 2) * dim + col + 1].red);
        g = (src[(dim - 2) * dim + col - 1].green + src[(dim - 2) * dim + col].green + src[(dim - 2) * dim + col + 1].green);
        //add the value
        dst[(dim - 3) * dim + col].blue += b;
        dst[(dim - 3) * dim + col].red += r;
        dst[(dim - 3) * dim + col].green += g;
        //add the value
        dst[(dim - 2) * dim + col].blue += b;
        dst[(dim - 2) * dim + col].red += r;
        dst[(dim - 2) * dim + col].green += g;
    }
    //row = dim - 1
    for (col = 1; col < dim - 1; ++col) {
        dst[(dim - 2) * dim + col].blue += (src[(dim - 1) * dim + col - 1].blue + src[(dim - 1) * dim + col].blue + src[(dim - 1) * dim + col + 1].blue);
        dst[(dim - 2) * dim + col].red += (src[(dim - 1) * dim + col - 1].red + src[(dim - 1) * dim + col].red + src[(dim - 1) * dim + col + 1].red);
        dst[(dim - 2) * dim + col].green += (src[(dim - 1) * dim + col - 1].green + src[(dim - 1) * dim + col].green + src[(dim - 1) * dim + col + 1].green);
    }

    int ii, jj,min_intensity, max_intensity, min_row, max_row, min_col, max_col;
    pixel loop_pixel;
    for (int i = 1; i < dim - 1; i++){
        for(int j = 1; j < dim -1; j++){
            min_intensity = 1024;
            max_intensity = -1;
            for(ii = i-1; ii <= i+1; ii++){
                for(jj = j-1 ; jj <= j+1 ; jj++) {
                    // check if smaller than min or higher than max and update
                    loop_pixel = src[calcIndex(ii, jj, dim)];
                    if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) <= min_intensity) {
                        min_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
                        min_row = ii;
                        min_col = jj;
                    }
                    if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) > max_intensity) {
                        max_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
                        max_row = ii;
                        max_col = jj;
                    }
                }
            }
            // filter out min and max
            dst[i*dim + j].red -= (src[min_row*dim + min_col].red + src[max_row*dim + max_col].red);
            dst[i*dim + j].blue -= (src[min_row*dim + min_col].blue + src[max_row*dim + max_col].blue);
            dst[i*dim + j].green -= (src[min_row*dim + min_col].green + src[max_row*dim + max_col].green);
            dst[i*dim + j].red /= 7;
            dst[i*dim + j].blue /= 7;
            dst[i*dim + j].green /= 7;
            if (dst[i*dim + j].blue < 0)
                dst[i*dim + j].blue = 0;
            else if (dst[i*dim + j].blue > 255)
                dst[i*dim + j].blue = 255;

            if (dst[i*dim + j].red < 0)
                dst[i*dim + j].red = 0;
            else if (dst[i*dim + j].red > 255)
                dst[i*dim + j].red = 255;

            if (dst[i*dim + j].green < 0)
                dst[i*dim + j].green = 0;
            else if (dst[i*dim + j].green > 255)
                dst[i*dim + j].green = 255;

        }
    }
}

void sharpKernelSmoothWithoutFilter(int dim, pixel* src, pixel* dst){
    int row, col;
    short b, r, g;

    //row 0
    for(col = 1; col < dim - 1; ++col){
        //initialize also
        dst[dim + col].blue = -1*(src[col - 1].blue + src[col].blue + src[col + 1].blue);
        dst[dim + col].red = -1*(src[col - 1].red + src[col].red + src[col + 1].red);
        dst[dim + col].green = -1*(src[col - 1].green + src[col].green + src[col + 1].green);
    }

    //row 1
    for(col = 1; col < dim - 1; ++col){
        b = src[dim + col -1].blue + src[dim + col].blue + src[dim + col + 1].blue;
        r = src[dim + col -1].red + src[dim + col].red + src[dim + col + 1].red;
        g = src[dim + col -1].green + src[dim + col].green + src[dim + col + 1].green;
        //add the value
        dst[dim + col].blue += -b + 10 * src[dim + col].blue;
        dst[dim + col].red += -r + 10 * src[dim + col].red;
        dst[dim + col].green += -g + 10 * src[dim + col].green;
        //initialize also
        dst[2*dim + col].blue = -b;
        dst[2*dim + col].red = -r;
        dst[2*dim + col].green = -g;
    }

    //rows 2 - (dim-3)
    for(row = 2; row <= dim - 3; ++row){
        for(col = 1; col < dim - 1; ++col){
            b = src[row*dim + col - 1].blue + src[row*dim + col].blue + src[row*dim + col + 1].blue;
            g = src[row*dim + col - 1].green + src[row*dim + col].green + src[row*dim + col + 1].green;
            r = src[row*dim + col - 1].red + src[row*dim + col].red + src[row*dim + col + 1].red;
            //add the value
            dst[(row-1)*dim + col].blue -= b;
            dst[(row-1)*dim + col].red -= r;
            dst[(row-1)*dim + col].green -= g;
            //check if the value between 0 and 255
            if (dst[(row-1)*dim + col].blue < 0)
                dst[(row-1)*dim + col].blue = 0;
            else if (dst[(row-1)*dim + col].blue > 255)
                dst[(row-1)*dim + col].blue = 255;

            if (dst[(row-1)*dim + col].red < 0)
                dst[(row-1)*dim + col].red = 0;
            else if (dst[(row-1)*dim + col].red > 255)
                dst[(row-1)*dim + col].red = 255;

            if (dst[(row-1)*dim + col].green < 0)
                dst[(row-1)*dim + col].green = 0;
            else if (dst[(row-1)*dim + col].green > 255)
                dst[(row-1)*dim + col].green = 255;
            //add the value + 10* the same
            dst[row*dim + col].blue += -b + 10*src[row*dim + col].blue;
            dst[row*dim + col].red += -r + 10*src[row*dim + col].red;
            dst[row*dim + col].green += -g + 10*src[row*dim + col].green;
            //initialize also
            dst[(row+1)*dim + col].blue = -b;
            dst[(row+1)*dim + col].red = -r;
            dst[(row+1)*dim + col].green = -g;
        }
    }
    //row = dim -2;
    for(col = 1; col < dim -1; ++col){
        b = src[(dim - 2)*dim + col -1].blue + src[(dim - 2)*dim + col].blue + src[(dim - 2)*dim + col + 1].blue;
        r = src[(dim - 2)*dim + col -1].red + src[(dim - 2)*dim + col].red + src[(dim - 2)*dim + col + 1].red;
        g = src[(dim - 2)*dim + col -1].green + src[(dim - 2)*dim + col].green + src[(dim - 2)*dim + col + 1].green;
        //add the value
        dst[(dim - 3)*dim + col].blue -= b;
        dst[(dim - 3)*dim + col].red -= r;
        dst[(dim - 3)*dim + col].green -= g;
        //close
        if (dst[(dim -3)*dim + col].blue < 0)
            dst[(dim -3)*dim + col].blue = 0;
        else if (dst[(dim -3)*dim + col].blue > 255)
            dst[(dim -3)*dim + col].blue = 255;

        if (dst[(dim -3)*dim + col].red < 0)
            dst[(dim -3)*dim + col].red = 0;
        else if (dst[(dim -3)*dim + col].red > 255)
            dst[(dim -3)*dim + col].red = 255;

        if (dst[(dim -3)*dim + col].green < 0)
            dst[(dim -3)*dim + col].green = 0;
        else if (dst[(dim -3)*dim + col].green > 255)
            dst[(dim -3)*dim + col].green = 255;
        //add the value
        dst[(dim - 2)*dim + col].blue += -b + 10 * src[(dim - 2)*dim + col].blue;
        dst[(dim - 2)*dim + col].red += -r + 10 * src[(dim - 2)*dim + col].red;
        dst[(dim - 2)*dim + col].green += -g + 10 * src[(dim - 2)*dim + col].green;
    }
    //row = dim - 1
    for(col = 1; col < dim -1; ++col){
        dst[(dim - 2)*dim + col].blue -= src[(dim - 1)*dim + col - 1].blue + src[(dim - 1)*dim + col].blue + src[(dim - 1)*dim + col + 1].blue;
        dst[(dim - 2)*dim + col].red -= src[(dim - 1)*dim + col - 1].red + src[(dim - 1)*dim + col].red + src[(dim - 1)*dim + col + 1].red;
        dst[(dim - 2)*dim + col].green -= src[(dim - 1)*dim + col - 1].green + src[(dim - 1)*dim + col].green + src[(dim - 1)*dim + col + 1].green;
        //close
        if (dst[(dim -2)*dim + col].blue < 0)
            dst[(dim -2)*dim + col].blue = 0;
        else if (dst[(dim -2)*dim + col].blue > 255)
            dst[(dim -2)*dim + col].blue = 255;

        if (dst[(dim -2)*dim + col].red < 0)
            dst[(dim -2)*dim + col].red = 0;
        else if (dst[(dim -2)*dim + col].red > 255)
            dst[(dim -2)*dim + col].red = 255;

        if (dst[(dim -2)*dim + col].green < 0)
            dst[(dim -2)*dim + col].green = 0;
        else if (dst[(dim -2)*dim + col].green > 255)
            dst[(dim -2)*dim + col].green = 255;

    }
}

void charsToPixels(Image *charsImg, pixel* pixels) {

	int row, col;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {

			pixels[row*n + col].red = (unsigned char) image->data[3*row*n + 3*col];
			pixels[row*n + col].green = (unsigned char) image->data[3*row*n + 3*col + 1];
			pixels[row*n + col].blue = (unsigned char) image->data[3*row*n + 3*col + 2];
		}
	}
}

void pixelsToChars(pixel* pixels, Image *charsImg) {

	int row, col;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {

			image->data[3*row*n + 3*col] = pixels[row*n + col].red;
			image->data[3*row*n + 3*col + 1] = pixels[row*n + col].green;
			image->data[3*row*n + 3*col + 2] = pixels[row*n + col].blue;
		}
	}
}

void copyPixels(pixel* src, pixel* dst) {

	int row, col;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {

			dst[row*n + col].red = src[row*n + col].red;
			dst[row*n + col].green = src[row*n + col].green;
			dst[row*n + col].blue = src[row*n + col].blue;
		}
	}
}

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

    pixel* pixelsImg = (pixel*)malloc(m*n*sizeof(pixel));
    pixel* backupOrg = (pixel*)malloc(m*n*sizeof(pixel));

    charsToPixels(image, pixelsImg);
    copyPixels(pixelsImg, backupOrg);

    smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

    pixelsToChars(pixelsImg, image);

    free(pixelsImg);
    free(backupOrg);
}

pixel* doConvolutionBlur9(Image *image) {

	pixel* pixelsImg = (pixel*)malloc(m*n*sizeof(pixel));
	pixel* backupOrg = (pixel*)malloc(m*n*sizeof(pixel));

	charsToPixels(image, pixelsImg);
	copyPixels(pixelsImg, backupOrg);

    blurKernelSmoothWithoutFilter9(m, backupOrg, pixelsImg);
	pixelsToChars(pixelsImg, image);

	free(backupOrg);
    return pixelsImg;
}

pixel* doConvolutionBlur7(Image *image) {

    pixel* pixelsImg = (pixel*)malloc(m*n*sizeof(pixel));
    pixel* backupOrg = (pixel*)malloc(m*n*sizeof(pixel));

    charsToPixels(image, pixelsImg);
    copyPixels(pixelsImg, backupOrg);

    blurKernelSmoothWithFilter7(m, backupOrg, pixelsImg);
    pixelsToChars(pixelsImg, image);

    //free(pixelsImg);
    free(backupOrg);
    return pixelsImg;
}

void doConvolutionSharp(Image *image, pixel* pixelsImg) {

    //pixel* pixelsImg = (pixel*)malloc(m*n*sizeof(pixel));
    pixel* backupOrg = (pixel*)malloc(m*n*sizeof(pixel));

    //charsToPixels(image, pixelsImg);
    copyPixels(pixelsImg, backupOrg);

    sharpKernelSmoothWithoutFilter(m, backupOrg, pixelsImg);
    pixelsToChars(pixelsImg, image);

    free(pixelsImg);
    free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

	/*
	* [1, 1, 1]
	* [1, 1, 1]
	* [1, 1, 1]
	*/
	//int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

	/*
	* [-1, -1, -1]
	* [-1, 9, -1]
	* [-1, -1, -1]
	*/
	//int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};
    pixel* pixelsImg;
	if (flag == '1') {	
		// blur image
        pixelsImg = doConvolutionBlur9(image);

		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);	

		// sharpen the resulting image
        doConvolutionSharp(image, pixelsImg);

		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);	
	} else {
		// apply extermum filtered kernel to blur image
        pixelsImg = doConvolutionBlur7(image);

		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image
        doConvolutionSharp(image, pixelsImg);

		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);	
	}
}

