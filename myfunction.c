#include <stdbool.h>
#include "writeBMP.h"
#include "readBMP.h"
#include <stdlib.h>
#define KERNEL_SIZE 3


typedef struct {
    short red;
    short green;
    short blue;
    short sum;
} pixel;

typedef struct {
    int min_loc;
    short min_val;
    int max_loc;
    short max_val;
} minMax;

int calcIndex(int i, int j, int n) {
	return ((i)*(n)+(j));
}

void blurKernelSmoothWithoutFilter9(int dim, pixel* src, pixel* dst){
    int row, col;
    int b, r, g;

    //row 0
    dst[0] = src[0];
    for(col = 1; col < dim - 1; ++col){
        dst[col] = src[col];
        //initialize also
        dst[dim + col].blue = (src[col - 1].blue + src[col].blue + src[col + 1].blue);
        dst[dim + col].red = (src[col - 1].red + src[col].red + src[col + 1].red);
        dst[dim + col].green = (src[col - 1].green + src[col].green + src[col + 1].green);
    }
    dst[dim -1] = src[dim -1];
    //row 1
    dst[dim] = src[dim];
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
    dst[2*dim -1] = src[2*dim -1];

    //rows 2 - (dim-3)
    for(row = 2; row <= dim - 3; ++row){
        dst[dim * row] = src[dim * row];
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
        dst[dim*row + dim -1] = src[dim*row + dim -1];
    }
    //row = dim -2;
    dst[(dim - 2)*dim] = src[(dim - 2)*dim];
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
    dst[(dim - 2)*dim + dim -1] = src[(dim - 2)*dim + dim - 1];
    //row = dim - 1
    dst[(dim - 1)*dim] = src[(dim - 1)*dim];
    for(col = 1; col < dim -1; ++col){
        dst[(dim - 1)*dim + col] = src[(dim - 1)*dim + col];
        dst[(dim - 2)*dim + col].blue += (src[(dim - 1)*dim + col - 1].blue + src[(dim - 1)*dim + col].blue + src[(dim - 1)*dim + col + 1].blue);
        dst[(dim - 2)*dim + col].red += (src[(dim - 1)*dim + col - 1].red + src[(dim - 1)*dim + col].red + src[(dim - 1)*dim + col + 1].red);
        dst[(dim - 2)*dim + col].green += (src[(dim - 1)*dim + col - 1].green + src[(dim - 1)*dim + col].green + src[(dim - 1)*dim + col + 1].green);
        //close
        dst[(dim - 2)*dim + col].blue /= 9;
        dst[(dim - 2)*dim + col].red /= 9;
        dst[(dim - 2)*dim + col].green /= 9;
        //copy the last row
    }
    dst[dim*dim - 1] = src[dim*dim -1];
}

void blurKernelSmoothWithFilter77(int dim, pixel* src, pixel* dst) {
    int row, col;
    int b, r, g;
    short a1, a2, a3;
    int cen;
    minMax* mM = (minMax *)malloc(dim*dim*sizeof(minMax));

    dst[0] = src[0];
    for (col = 1; col < dim - 1; ++col) {
        //initialize also
        dst[col] = src[col];
        dst[dim + col].blue = (src[col - 1].blue + src[col].blue + src[col + 1].blue);
        dst[dim + col].red = (src[col - 1].red + src[col].red + src[col + 1].red);
        dst[dim + col].green = (src[col - 1].green + src[col].green + src[col + 1].green);

        //sort sum

        cen = col;
        a1 = src[cen -1].sum;
        a2 = src[cen].sum;
        a3 = src[cen + 1].sum;

        if (a1 >= a2){
            if (a1 >= a3){
                mM[cen].max_loc = cen - 1;
                mM[cen].max_val = a1;
                if (a2 >= a3) {
                    mM[cen].min_loc = cen + 1;
                    mM[cen].min_val = a3;
                }
                else{
                    mM[cen].min_loc = cen;
                    mM[cen].min_val = a2;
                }
            }
            else{
                mM[cen].max_loc = cen + 1;
                mM[cen].max_val = a3;
                mM[cen].min_loc = cen;
                mM[cen].min_val = a2;
            }
        }
        else{
            if (a2 >= a3){
                mM[cen].max_loc = cen;
                mM[cen].max_val = a2;
                if (a3 > a1) {
                    mM[cen].min_loc = cen - 1;
                    mM[cen].min_val = a1;
                }
                else {
                    mM[cen].min_loc = cen + 1;
                    mM[cen].min_val = a3;
                }
            }
            else{
                mM[cen].max_loc = cen + 1;
                mM[cen].max_val = a3;
                mM[cen].min_loc = cen - 1;
                mM[cen].min_val = a1;
            }
        }
    }
    dst[dim - 1] = src[dim -1];
    //row 1
    dst[dim] = src[dim];
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

        cen = dim + col;
        a1 = src[cen -1].sum;
        a2 = src[cen].sum;
        a3 = src[cen + 1].sum;

        if (a1 >= a2){
            if (a1 >= a3){
                mM[cen].max_loc = cen - 1;
                mM[cen].max_val = a1;
                if (a2 >= a3) {
                    mM[cen].min_loc = cen + 1;
                    mM[cen].min_val = a3;
                }
                else{
                    mM[cen].min_loc = cen;
                    mM[cen].min_val = a2;
                }
            }
            else{
                mM[cen].max_loc = cen + 1;
                mM[cen].max_val = a3;
                mM[cen].min_loc = cen;
                mM[cen].min_val = a2;
            }
        }
        else{
            if (a2 >= a3){
                mM[cen].max_loc = cen;
                mM[cen].max_val = a2;
                if (a3 > a1) {
                    mM[cen].min_loc = cen - 1;
                    mM[cen].min_val = a1;
                }
                else {
                    mM[cen].min_loc = cen + 1;
                    mM[cen].min_val = a3;
                }
            }
            else{
                mM[cen].max_loc = cen + 1;
                mM[cen].max_val = a3;
                mM[cen].min_loc = cen - 1;
                mM[cen].min_val = a1;
            }
        }

    }
    dst[2*dim - 1] = src[2*dim -1];
    //rows 2 - (dim-3)
    for (row = 2; row <= dim - 3; ++row) {
        dst[row*dim] = src[row*dim];
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

            cen = row*dim + col;
            a1 = src[cen -1].sum;
            a2 = src[cen].sum;
            a3 = src[cen + 1].sum;

            if (a1 >= a2){
                if (a1 >= a3){
                    mM[cen].max_loc = cen - 1;
                    mM[cen].max_val = a1;
                    if (a2 >= a3) {
                        mM[cen].min_loc = cen + 1;
                        mM[cen].min_val = a3;
                    }
                    else{
                        mM[cen].min_loc = cen;
                        mM[cen].min_val = a2;
                    }
                }
                else{
                    mM[cen].max_loc = cen + 1;
                    mM[cen].max_val = a3;
                    mM[cen].min_loc = cen;
                    mM[cen].min_val = a2;
                }
            }
            else{
                if (a2 >= a3){
                    mM[cen].max_loc = cen;
                    mM[cen].max_val = a2;
                    if (a3 > a1) {
                        mM[cen].min_loc = cen - 1;
                        mM[cen].min_val = a1;
                    }
                    else {
                        mM[cen].min_loc = cen + 1;
                        mM[cen].min_val = a3;
                    }
                }
                else{
                    mM[cen].max_loc = cen + 1;
                    mM[cen].max_val = a3;
                    mM[cen].min_loc = cen - 1;
                    mM[cen].min_val = a1;
                }
            }
        }
        dst[row*dim + dim -1] = src[row*dim + dim - 1];
    }
    //row = dim -2;
    dst[(dim - 2) * dim] = src[(dim - 2) * dim];
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

        cen = (dim - 2) * dim + col;
        a1 = src[cen -1].sum;
        a2 = src[cen].sum;
        a3 = src[cen + 1].sum;

        if (a1 >= a2){
            if (a1 >= a3){
                mM[cen].max_loc = cen - 1;
                mM[cen].max_val = a1;
                if (a2 >= a3) {
                    mM[cen].min_loc = cen + 1;
                    mM[cen].min_val = a3;
                }
                else{
                    mM[cen].min_loc = cen;
                    mM[cen].min_val = a2;
                }
            }
            else{
                mM[cen].max_loc = cen + 1;
                mM[cen].max_val = a3;
                mM[cen].min_loc = cen;
                mM[cen].min_val = a2;
            }
        }
        else{
            if (a2 >= a3){
                mM[cen].max_loc = cen;
                mM[cen].max_val = a2;
                if (a3 > a1) {
                    mM[cen].min_loc = cen - 1;
                    mM[cen].min_val = a1;
                }
                else {
                    mM[cen].min_loc = cen + 1;
                    mM[cen].min_val = a3;
                }
            }
            else{
                mM[cen].max_loc = cen + 1;
                mM[cen].max_val = a3;
                mM[cen].min_loc = cen - 1;
                mM[cen].min_val = a1;
            }
        }
    }
    dst[(dim - 2) * dim + dim - 1] = src[(dim - 2) * dim + dim -1];
    //row = dim - 1
    dst[(dim -1)*dim] = src[(dim - 1)*dim];
    for (col = 1; col < dim - 1; ++col) {
        dst[(dim -1)*dim + col] = src[(dim -1)*dim + col];
        dst[(dim - 2) * dim + col].blue += (src[(dim - 1) * dim + col - 1].blue + src[(dim - 1) * dim + col].blue + src[(dim - 1) * dim + col + 1].blue);
        dst[(dim - 2) * dim + col].red += (src[(dim - 1) * dim + col - 1].red + src[(dim - 1) * dim + col].red + src[(dim - 1) * dim + col + 1].red);
        dst[(dim - 2) * dim + col].green += (src[(dim - 1) * dim + col - 1].green + src[(dim - 1) * dim + col].green + src[(dim - 1) * dim + col + 1].green);

        cen = (dim - 1) * dim + col;
        a1 = src[cen -1].sum;
        a2 = src[cen].sum;
        a3 = src[cen + 1].sum;

        if (a1 >= a2){
            if (a1 >= a3){
                mM[cen].max_loc = cen - 1;
                mM[cen].max_val = a1;
                if (a2 >= a3) {
                    mM[cen].min_loc = cen + 1;
                    mM[cen].min_val = a3;
                }
                else{
                    mM[cen].min_loc = cen;
                    mM[cen].min_val = a2;
                }
            }
            else{
                mM[cen].max_loc = cen + 1;
                mM[cen].max_val = a3;
                mM[cen].min_loc = cen;
                mM[cen].min_val = a2;
            }
        }
        else{
            if (a2 >= a3){
                mM[cen].max_loc = cen;
                mM[cen].max_val = a2;
                if (a3 > a1) {
                    mM[cen].min_loc = cen - 1;
                    mM[cen].min_val = a1;
                }
                else {
                    mM[cen].min_loc = cen + 1;
                    mM[cen].min_val = a3;
                }
            }
            else{
                mM[cen].max_loc = cen + 1;
                mM[cen].max_val = a3;
                mM[cen].min_loc = cen - 1;
                mM[cen].min_val = a1;
            }
        }
    }
    dst[dim*dim -1 ] = src[dim * dim - 1];
    int minLocation, maxLocation;
    short m1, m2, m3, M1, M2, M3;
    for (int i = 1; i < dim - 1; i++){
        for(int j = 1; j < dim -1; j++){
            //initialize
            m1 = mM[(i - 1)*dim + j].min_val;
            m2 = mM[i*dim + j].min_val;
            m3 = mM[(i + 1)*dim + j].min_val;

            M1 = mM[(i - 1)*dim + j].max_val;
            M2 = mM[i*dim + j].max_val;
            M3 = mM[(i + 1)*dim + j].max_val;

            //find min
            if (m1 < m2){
                if (m1 < m3)
                    minLocation = mM[(i - 1)*dim + j].min_loc;
                else
                    minLocation = mM[(i + 1)*dim + j].min_loc;
            }
            else{
                if (m3 <= m2)
                    minLocation = mM[(i + 1)*dim + j].min_loc;
                else
                    minLocation = mM[i*dim + j].min_loc;
            }

            //find max
            if (M1 >= M2){
                if (M1 >= M3)
                    maxLocation = mM[(i - 1)*dim + j].max_loc;
                else
                    maxLocation = mM[(i + 1)*dim + j].max_loc;
            }
            else{
                if (M3 > M2)
                    maxLocation = mM[(i + 1)*dim + j].max_loc;
                else
                    maxLocation = mM[i*dim + j].max_loc;
            }

            // filter out min and max
            dst[i*dim + j].red -= (src[minLocation].red + src[maxLocation].red);
            dst[i*dim + j].blue -= (src[minLocation].blue + src[maxLocation].blue);
            dst[i*dim + j].green -= (src[minLocation].green + src[maxLocation].green);

            //norm
            dst[i*dim + j].red /= 7;
            dst[i*dim + j].blue /= 7;
            dst[i*dim + j].green /= 7;

            //check valid
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
    free(mM);
}

void sharpKernelSmoothWithoutFilter(int dim, pixel* src, pixel* dst){
    int row, col;
    short b, r, g;

    //row 0
    dst[0] = src[0];
    for(col = 1; col < dim - 1; ++col){
        //initialize also
        dst[col] = src[col];
        dst[dim + col].blue = -1*(src[col - 1].blue + src[col].blue + src[col + 1].blue);
        dst[dim + col].red = -1*(src[col - 1].red + src[col].red + src[col + 1].red);
        dst[dim + col].green = -1*(src[col - 1].green + src[col].green + src[col + 1].green);
    }
    dst[dim -1] = src[dim - 1];
    //row 1
    dst[dim] = src[dim];
    for(col = 1; col < dim - 1; ++col){
        b = src[dim + col -1].blue + src[dim + col].blue + src[dim + col + 1].blue;
        r = src[dim + col -1].red + src[dim + col].red + src[dim + col + 1].red;
        g = src[dim + col -1].green + src[dim + col].green + src[dim + col + 1].green;
        //add the value
        dst[dim + col].blue += 10 * src[dim + col].blue - b;
        dst[dim + col].red += 10 * src[dim + col].red - r;
        dst[dim + col].green += 10 * src[dim + col].green - g;
        //initialize also
        dst[2*dim + col].blue = -b;
        dst[2*dim + col].red = -r;
        dst[2*dim + col].green = -g;
    }
    dst[2*dim -1] = src[2*dim -1];
    short bb, rr, gg;
    //rows 2 - (dim-3)
    for(row = 2; row <= dim - 3; ++row){
        dst[row*dim] = src[row*dim];
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
            dst[row*dim + col].blue += 10*src[row*dim + col].blue - b;
            dst[row*dim + col].red += 10*src[row*dim + col].red - r;
            dst[row*dim + col].green += 10*src[row*dim + col].green - g;

            //initialize also
            dst[(row+1)*dim + col].blue = -b;
            dst[(row+1)*dim + col].red = -r;
            dst[(row+1)*dim + col].green = -g;
        }
        dst[row*dim + dim - 1] = src[row*dim + dim -1];
    }
    //row = dim -2;
    dst[(dim - 2)*dim] = src[(dim - 2)*dim];
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
    dst[(dim - 2)*dim + dim -1] = src[(dim - 2)*dim + dim - 1];
    //row = dim - 1
    dst[(dim -1)*dim] = src[(dim - 1)*dim];
    for(col = 1; col < dim -1; ++col){
        dst[(dim -1)*dim + col] = src[(dim - 1)*dim + col];
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
    dst[dim*dim - 1] = src[dim*dim - 1];
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

void charsToPixels7(Image *charsImg, pixel* pixels) {

    int row, col;
    for (row = 0 ; row < m ; row++) {
        for (col = 0 ; col < n ; col++) {

            pixels[row*n + col].red = (unsigned char) image->data[3*row*n + 3*col];
            pixels[row*n + col].green = (unsigned char) image->data[3*row*n + 3*col + 1];
            pixels[row*n + col].blue = (unsigned char) image->data[3*row*n + 3*col + 2];
            pixels[row*n + col].sum = pixels[row*n + col].blue + pixels[row*n + col].green + pixels[row*n + col].red;
        }
    }
}

void pixelsToChars(pixel* pixels, Image *charsImg) {

	int row, col;
    int rowNcol = 0, row3Ncol = 0;
	for (row = 0 ; row < m ; row++) {
		for (col = 0 ; col < n ; col++) {

			image->data[3*row*n + 3*col] = pixels[row*n + col].red;
			image->data[3*row*n + 3*col + 1] = pixels[row*n + col].green;
			image->data[3*row*n + 3*col + 2] = pixels[row*n + col].blue;
		}
	}
}

pixel* doConvolutionBlur9(Image *image) {

	pixel* pixelsImg = (pixel*)malloc(m*n*sizeof(pixel));
	pixel* backupOrg = (pixel*)malloc(m*n*sizeof(pixel));

	charsToPixels(image, backupOrg);
	//copyPixels(pixelsImg, backupOrg);

    blurKernelSmoothWithoutFilter9(m, backupOrg, pixelsImg);
	pixelsToChars(pixelsImg, image);

	free(backupOrg);
    return pixelsImg;
}

pixel* doConvolutionBlur7(Image *image) {

    pixel* pixelsImg = (pixel*)malloc(m*n*sizeof(pixel));
    pixel* backupOrg = (pixel*)malloc(m*n*sizeof(pixel));

    charsToPixels7(image, backupOrg);
    //copyPixels(pixelsImg, backupOrg);

    blurKernelSmoothWithFilter77(m, backupOrg, pixelsImg);
    pixelsToChars(pixelsImg, image);

    //free(pixelsImg);
    free(backupOrg);
    return pixelsImg;
}

void doConvolutionSharp(Image *image, pixel* pixelsImg) {

    pixel* backupOrg = (pixel*)malloc(m*n*sizeof(pixel));

    //charsToPixels(image, pixelsImg);
    //copyPixels(pixelsImg, backupOrg);

    sharpKernelSmoothWithoutFilter(m,pixelsImg, backupOrg);
    pixelsToChars(backupOrg, image);

    free(pixelsImg);
    free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

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

