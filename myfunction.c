//Sahar Rofe 209275114
#include <stdbool.h>

/*
 * this struct is the pixel!!!
 * I use the sum only with the filter.
 */
typedef struct {
    short red;
    short green;
    short blue;
    short sum;
} pixel;

/*
 * I use this struct only in the filter
 */
typedef struct {
    int min_loc;
    int max_loc;
    short min_val;
    short max_val;
} minMax;


/*
 * this function replace the smooth with the regular blur without filter
 */
void blurKernelSmoothWithoutFilter9(int dim, pixel* src, pixel* dst){
    int row, col;
    int b, r, g;
    int dimRowCol, dimRowColM1, dimRowColP1;
    register unsigned long d;
    register unsigned long s1, s2, s3;
    register unsigned long ss1, ss2, ss3;
    int mod = (dim-1) % 4;

    //row 0
    dst[0] = src[0];
    for(col = 1; col < dim - 1; ++col){
        *(unsigned long*)(dst + col) = *(unsigned long*)(src + col);
        //initialize also
        s1 = *(unsigned long*)(src + col - 1);
        s2 = *(unsigned long*)(src + col);
        s3 = *(unsigned long*)(src + col + 1);
        r = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
        s1 >>= 16;
        s2 >>= 16;
        s3 >>= 16;
        g = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
        s1 >>= 16;
        s2 >>= 16;
        s3 >>= 16;
        b = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
        dst[dim + col].blue = b;
        dst[dim + col].red = r;
        dst[dim + col].green = g;

    }
    dst[dim -1] = src[dim -1];
    //row 1
    dst[dim] = src[dim];
    dimRowColP1 = 2*dim;
    dimRowCol = dim;
    for(col = 1; col < dim - 1; ++col){
        ++dimRowCol;
        ++dimRowColP1;
        s1 = *(unsigned long*)(src + dimRowCol - 1);
        s2 = *(unsigned long*)(src + dimRowCol);
        s3 = *(unsigned long*)(src + dimRowCol + 1);
        r = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
        s1 >>= 16;
        s2 >>= 16;
        s3 >>= 16;
        g = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
        s1 >>= 16;
        s2 >>= 16;
        s3 >>= 16;
        b = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
        d = ((((long)(unsigned short)(b)) << 16) << 16) +
            ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
        *((long *)(dst + dimRowCol)) += d;
        *((long *)(dst + dimRowColP1)) = d;
    }
    ++dimRowCol;
    dst[dimRowCol] = src[dimRowCol];

    //rows 2 - (dim-3)
    dimRowColM1 = dim;
    dimRowCol = 2*dim;
    dimRowColP1 = 3*dim;
    for(row = 2; row <= dim - 3; ++row){
        *((unsigned long*)(dst + dimRowCol)) = *((unsigned long*)(src + dimRowCol));
        for(col = 1; col < mod; ++col){
            ++dimRowCol;
            ++dimRowColP1;
            ++dimRowColM1;
            s1 = *(unsigned long*)(src + dimRowCol - 1);
            s2 = *(unsigned long*)(src + dimRowCol);
            s3 = *(unsigned long*)(src + dimRowCol + 1);
            r = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            g = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            b = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;

            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColM1)) += d;
            *((long *)(dst + dimRowCol)) += d;
            *((long *)(dst + dimRowColP1)) = d;

            dst[dimRowColM1].blue /= 9;
            dst[dimRowColM1].red /= 9;
            dst[dimRowColM1].green /= 9;
        }

        for(; col < dim - 1; col += 4){
          ++dimRowCol;
          ++dimRowColP1;
          ++dimRowColM1;
          ss1 = *(unsigned long*)(src + dimRowCol - 1);
          ss2 = *(unsigned long*)(src + dimRowCol);
          ss3 = *(unsigned long*)(src + dimRowCol + 1);
          s1 = ss1;
          s2 = ss2;
          s3 = ss3;
          r = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
          s1 >>= 16;
          s2 >>= 16;
          s3 >>= 16;
          g = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
          s1 >>= 16;
          s2 >>= 16;
          s3 >>= 16;
          b = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;

          d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
          *((long *)(dst + dimRowColM1)) += d;
          *((long *)(dst + dimRowCol)) += d;
          *((long *)(dst + dimRowColP1)) = d;

          dst[dimRowColM1].blue /= 9;
          dst[dimRowColM1].red /= 9;
          dst[dimRowColM1].green /= 9;
          ////////////////////////////////////////////
            ++dimRowCol;
            ++dimRowColP1;
            ++dimRowColM1;
            ss1 = ss2;
            ss2 = ss3;
            ss3 = *(unsigned long*)(src + dimRowCol + 1);
            s1 = ss1;
            s2 = ss2;
            s3 = ss3;
            r = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            g = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            b = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;

            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColM1)) += d;
            *((long *)(dst + dimRowCol)) += d;
            *((long *)(dst + dimRowColP1)) = d;

            dst[dimRowColM1].blue /= 9;
            dst[dimRowColM1].red /= 9;
            dst[dimRowColM1].green /= 9;
            ///////////////////////////////////////
            ++dimRowCol;
            ++dimRowColP1;
            ++dimRowColM1;
            ss1 = ss2;
            ss2 = ss3;
            ss3 = *(unsigned long*)(src + dimRowCol + 1);
            s1 = ss1;
            s2 = ss2;
            s3 = ss3;
            r = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            g = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            b = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;

            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColM1)) += d;
            *((long *)(dst + dimRowCol)) += d;
            *((long *)(dst + dimRowColP1)) = d;

            dst[dimRowColM1].blue /= 9;
            dst[dimRowColM1].red /= 9;
            dst[dimRowColM1].green /= 9;
            //////////////////////////////////////////
            ++dimRowCol;
            ++dimRowColP1;
            ++dimRowColM1;
            ss1 = ss2;
            ss2 = ss3;
            s3 = *(unsigned long*)(src + dimRowCol + 1);
            s1 = ss1;
            s2 = ss2;
            r = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            g = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            b = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;

            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColM1)) += d;
            *((long *)(dst + dimRowCol)) += d;
            *((long *)(dst + dimRowColP1)) = d;

            dst[dimRowColM1].blue /= 9;
            dst[dimRowColM1].red /= 9;
            dst[dimRowColM1].green /= 9;
        }
        dst[dimRowCol + 1] = src[dimRowCol + 1];
        dimRowCol += 2;
        dimRowColM1 += 2;
        dimRowColP1 += 2;
    }
    //row = dim -2;
    //dimRow = dim*(dim -2);
    dst[dimRowCol] = src[dimRowCol];
    for(col = 1; col < dim -1; ++col){
        ++dimRowCol;
        ++dimRowColM1;
        b = (src[dimRowCol -1].blue + src[dimRowCol].blue + src[dimRowCol + 1].blue);
        r = (src[dimRowCol -1].red + src[dimRowCol].red + src[dimRowCol + 1].red);
        g = (src[dimRowCol -1].green + src[dimRowCol].green + src[dimRowCol + 1].green);
        //add the value
        dst[dimRowColM1].blue += b;
        dst[dimRowColM1].red += r;
        dst[dimRowColM1].green += g;
        //close
        dst[dimRowColM1].blue /= 9;
        dst[dimRowColM1].red /= 9;
        dst[dimRowColM1].green /= 9;
        //add the value
        dst[dimRowCol].blue += b;
        dst[dimRowCol].red += r;
        dst[dimRowCol].green += g;
    }
    ++dimRowCol;
    dst[dimRowCol] = src[dimRowCol];
    //row = dim - 1
    ++dimRowCol;
    dimRowColM1 += 2;
    dst[dimRowCol] = src[dimRowCol];
    for(col = 1; col < dim -1; ++col){
        ++dimRowCol;
        ++dimRowColM1;
        dst[dimRowCol] = src[dimRowCol];
        dst[dimRowColM1].blue += (src[dimRowCol - 1].blue + src[dimRowCol].blue + src[dimRowCol + 1].blue);
        dst[dimRowColM1].red += (src[dimRowCol - 1].red + src[dimRowCol].red + src[dimRowCol + 1].red);
        dst[dimRowColM1].green += (src[dimRowCol - 1].green + src[dimRowCol].green + src[dimRowCol + 1].green);
        //close
        dst[dimRowColM1].blue /= 9;
        dst[dimRowColM1].red /= 9;
        dst[dimRowColM1].green /= 9;
        //copy the last row
    }
    ++dimRowCol;
    dst[dimRowCol] = src[dimRowCol];
}

//this function replace the smooth with the spacial blur with filter
void blurKernelSmoothWithFilter77(int dim, pixel* src, pixel* dst) {
    int row, col;
    int b, r, g;
    short a1, a2, a3;
    int cen;
    int dimRowCol, dimRowColM1, dimRowColP1;
    short dstR, dstB, dstG;
    //in this array I calculate the min and max value and location
    minMax* mM = (minMax *)malloc(dim*dim*sizeof(minMax));
    register unsigned long s1;
    register unsigned long s2;
    register unsigned long s3;
    register unsigned long d;

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
    dimRowCol = dim;
    dimRowColP1 = 2*dim;
    for (col = 1; col < dim - 1; ++col) {
        ++dimRowCol;
        ++dimRowColP1;
        b = (src[dimRowCol - 1].blue + src[dimRowCol].blue + src[dimRowCol + 1].blue);
        r = (src[dimRowCol - 1].red + src[dimRowCol].red + src[dimRowCol + 1].red);
        g = (src[dimRowCol - 1].green + src[dimRowCol].green + src[dimRowCol + 1].green);
        //add the value
        dst[dimRowCol].blue += b;
        dst[dimRowCol].red += r;
        dst[dimRowCol].green += g;
        //initialize also
        dst[dimRowColP1].blue = b;
        dst[dimRowColP1].red = r;
        dst[dimRowColP1].green = g;

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
    dimRowColM1 = dim;
    dimRowCol = 2*dim;
    dimRowColP1 = 3*dim;
    for (row = 2; row <= dim - 3; ++row) {
        dst[dimRowCol] = src[dimRowCol];
        for (col = 1; col < dim - 1; ++col) {
            ++dimRowColM1;
            ++dimRowCol;
            ++dimRowColP1;
            s1 = *(unsigned long*)(src + dimRowCol - 1);
            s2 = *(unsigned long*)(src + dimRowCol);
            s3 = *(unsigned long*)(src + dimRowCol + 1);
            r = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            g = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            b = (unsigned char)s1 + (unsigned char)s2 + (unsigned char)s3;
            s1 >>= 16;
            s2 >>= 16;
            s3 >>= 16;
            a1 = (unsigned short)s1;
            a2 = (unsigned short)s2;
            a3 = (unsigned short)s3;
            //initialize also
            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColM1)) += d;
            *((long *)(dst + dimRowCol)) += d;
            *((long *)(dst + dimRowColP1)) = d;

            if (a1 >= a2){
                if (a1 >= a3){
                    mM[dimRowCol].max_loc = dimRowCol - 1;
                    mM[dimRowCol].max_val = a1;
                    if (a2 >= a3) {
                        mM[dimRowCol].min_loc = dimRowCol + 1;
                        mM[dimRowCol].min_val = a3;
                    }
                    else{
                        mM[dimRowCol].min_loc = dimRowCol;
                        mM[dimRowCol].min_val = a2;
                    }
                }
                else{
                    mM[dimRowCol].max_loc = dimRowCol + 1;
                    mM[dimRowCol].max_val = a3;
                    mM[dimRowCol].min_loc = dimRowCol;
                    mM[dimRowCol].min_val = a2;
                }
            }
            else{
                if (a2 >= a3){
                    mM[dimRowCol].max_loc = dimRowCol;
                    mM[dimRowCol].max_val = a2;
                    if (a3 > a1) {
                        mM[dimRowCol].min_loc = dimRowCol - 1;
                        mM[dimRowCol].min_val = a1;
                    }
                    else {
                        mM[dimRowCol].min_loc = dimRowCol + 1;
                        mM[dimRowCol].min_val = a3;
                    }
                }
                else{
                    mM[dimRowCol].max_loc = dimRowCol + 1;
                    mM[dimRowCol].max_val = a3;
                    mM[dimRowCol].min_loc = dimRowCol - 1;
                    mM[dimRowCol].min_val = a1;
                }
            }
        }
        ++dimRowCol;
        dst[dimRowCol] = src[dimRowCol];
        ++dimRowCol;
        dimRowColP1 += 2;
        dimRowColM1 += 2;
    }
    //row = dim -2;
    dst[dimRowCol] = src[dimRowCol];
    for (col = 1; col < dim - 1; ++col) {
        ++dimRowCol;
        ++dimRowColP1;
        ++dimRowColM1;
        b = (src[dimRowCol - 1].blue + src[dimRowCol].blue + src[dimRowCol + 1].blue);
        r = (src[dimRowCol - 1].red + src[dimRowCol].red + src[dimRowCol + 1].red);
        g = (src[dimRowCol - 1].green + src[dimRowCol].green + src[dimRowCol + 1].green);
        //add the value
        dst[dimRowColM1].blue += b;
        dst[dimRowColM1].red += r;
        dst[dimRowColM1].green += g;
        //add the value
        dst[dimRowCol].blue += b;
        dst[dimRowCol].red += r;
        dst[dimRowCol].green += g;

        cen = dimRowCol;
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
    ++dimRowCol;
    dst[dimRowCol] = src[dimRowCol];
    //row = dim - 1
    ++dimRowCol;
    dimRowColM1 += 2;
    dst[dimRowCol] = src[dimRowCol];
    for (col = 1; col < dim - 1; ++col) {
        ++dimRowCol;
        ++dimRowColM1;
        dst[dimRowCol] = src[dimRowCol];
        dst[(dim - 2) * dim + col].blue += (src[dimRowCol - 1].blue + src[dimRowCol].blue + src[dimRowCol + 1].blue);
        dst[(dim - 2) * dim + col].red += (src[dimRowCol - 1].red + src[dimRowCol].red + src[dimRowCol + 1].red);
        dst[(dim - 2) * dim + col].green += (src[dimRowCol - 1].green + src[dimRowCol].green + src[dimRowCol + 1].green);

        cen = dimRowCol;
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
    int i,j;
    dimRowColM1 = 0;
    dimRowCol = dim;
    dimRowColP1 = 2*dim;
    for (i = 1; i < dim - 1; i++){
        for(j = 1; j < dim -1; j++){
            ++dimRowColM1;
            ++dimRowCol;
            ++dimRowColP1;
            m1 = mM[dimRowColM1].min_val;
            m2 = mM[dimRowCol].min_val;
            m3 = mM[dimRowColP1].min_val;

            M1 = mM[dimRowColM1].max_val;
            M2 = mM[dimRowCol].max_val;
            M3 = mM[dimRowColP1].max_val;

            //find min
            if (m1 < m2){
                if (m1 < m3)
                    minLocation = mM[dimRowColM1].min_loc;
                else
                    minLocation = mM[dimRowColP1].min_loc;
            }
            else{
                if (m3 <= m2)
                    minLocation = mM[dimRowColP1].min_loc;
                else
                    minLocation = mM[dimRowCol].min_loc;
            }

            //find max
            if (M1 >= M2){
                if (M1 >= M3)
                    maxLocation = mM[dimRowColM1].max_loc;
                else
                    maxLocation = mM[dimRowColP1].max_loc;
            }
            else{
                if (M3 > M2)
                    maxLocation = mM[dimRowColP1].max_loc;
                else
                    maxLocation = mM[dimRowCol].max_loc;
            }

            // filter out min and max
            dstR = dst[dimRowCol].red;
            dstB = dst[dimRowCol].blue;
            dstG = dst[dimRowCol].green;

            dstR -= (src[minLocation].red + src[maxLocation].red);
            dstB -= (src[minLocation].blue + src[maxLocation].blue);
            dstG -= (src[minLocation].green + src[maxLocation].green);

            //norm
            dstR /= 7;
            dstB /= 7;
            dstG /= 7;

            //check valid
            if (dstB < 0)
                dstB = 0;
            else if (dstB > 255)
                dstB = 255;

            if (dstR < 0)
                dstR = 0;
            else if (dstR > 255)
                dstR = 255;

            if (dstG < 0)
                dstG = 0;
            else if (dstG > 255)
                dstG = 255;

            dst[dimRowCol].red = dstR;
            dst[dimRowCol].blue = dstB;
            dst[dimRowCol].green = dstG;
        }
        dimRowCol += 2;
        dimRowColP1 += 2;
        dimRowColM1 += 2;
    }
    free(mM);
}


//this function replace the smooth with sharp
void sharpKernelSmoothWithoutFilter(int dim, pixel* src, pixel* dst){

    // init var
    int row, col;
    short b, r, g;
    short dstB, dstR, dstG;
    int dimRow, dimRowM1;
    int dimRowCol, dimRowColM1, dimRowColP1;
    short srcB2, srcR2, srcG2, srcB1, srcG1, srcR1;
    int mod = (dim - 1) % 4;
    register unsigned long s1;
    register unsigned long s2;
    register unsigned long s3;
    register unsigned long ss1;
    register unsigned long ss2;
    register unsigned long ss3;
    register unsigned long d;

    // row 0

    for(col = 1; col < mod ; ++col){
        dimRow = dim + col;
        dst[dimRow].blue = -1*(src[col - 1].blue + src[col].blue + src[col + 1].blue);
        dst[dimRow].red = -1*(src[col - 1].red + src[col].red + src[col + 1].red);
        dst[dimRow].green = -1*(src[col - 1].green + src[col].green + src[col + 1].green);
    }

    for(; col < dim - 1; col += 4){
        dimRow = dim + col;
        srcB1 = src[col + 1].blue;
        srcG1 = src[col + 1].green;
        srcR1 = src[col + 1].red;

        srcB2 = src[col + 2].blue;
        srcG2 = src[col + 2].green;
        srcR2 = src[col + 2].red;

        dst[dimRow].blue = -1*(src[col - 1].blue + src[col].blue + srcB1);
        dst[dimRow].red = -1*(src[col - 1].red + src[col].red + srcR1);
        dst[dimRow].green = -1*(src[col - 1].green + src[col].green + srcG1);

        dst[dimRow + 1].blue = -1*(src[col].blue + srcB1 + srcB2);
        dst[dimRow + 1].red = -1*(src[col].red + srcR1 + srcR2);
        dst[dimRow + 1].green = -1*(src[col].green + srcG1 + srcG2);

        dst[dimRow + 2].blue = -1*(srcB1 + srcB2 + src[col + 3].blue);
        dst[dimRow + 2].red = -1*(srcR1 + srcR2 + src[col + 3].red);
        dst[dimRow + 2].green = -1*(srcG1 + srcG2 + src[col + 3].green);

        dst[dimRow + 3].blue = -1*(srcB2 + src[col + 3].blue + src[col + 4].blue);
        dst[dimRow + 3].red = -1*(srcR2 + src[col + 3].red + src[col + 4].red);
        dst[dimRow + 3].green = -1*(srcG2 + src[col + 3].green + src[col + 4].green);
    }

    //row 1
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

    //rows 2 - (dim-3)
    dimRowColM1 = dim;
    dimRowCol = 2*dim;
    dimRowColP1 = 3*dim;
    for(row = 2; row <= dim - 3; ++row){
        for(col = 1; col < mod; col++){
            ++dimRowCol;
            ++dimRowColM1;
            ++dimRowColP1;
            s1 = *(unsigned long*)(src + dimRowCol - 1);
            s2 = *(unsigned long*)(src + dimRowCol);
            s3 = *(unsigned long*)(src + dimRowCol + 1);
            srcR1 = (unsigned char)s2;
            s2 >>= 16;
            srcG1 = (unsigned char)s2;
            s2 >>= 16;
            srcB1 = (unsigned char)s2;
            r = -(unsigned char)s1 - srcR1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            g = -(unsigned char)s1 - srcG1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            b = -(unsigned char)s1 - srcB1 - (unsigned char)s3;

            dst[dimRowColM1].blue += b;
            dst[dimRowColM1].red += r;
            dst[dimRowColM1].green += g;

            //check if the value between 0 and 255
            dstB = dst[dimRowColM1].blue;
            dstR = dst[dimRowColM1].red;
            dstG = dst[dimRowColM1].green;

            if (dstB < 0)
                dst[dimRowColM1].blue = 0;
            else if (dstB > 255)
                dst[dimRowColM1].blue = 255;
            if (dstR < 0)
                dst[dimRowColM1].red = 0;
            else if (dstR > 255)
                dst[dimRowColM1].red = 255;
            if (dstG < 0)
                dst[dimRowColM1].green = 0;
            else if (dstG > 255)
                dst[dimRowColM1].green = 255;

            //add the value + 10* the same
            dst[dimRowCol].blue += 10*srcB1 + b;
            dst[dimRowCol].red += 10*srcR1 + r;
            dst[dimRowCol].green += 10*srcG1 + g;

            //initialize also
            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColP1)) = d;
        }

        for(; col < dim - 1; col += 4){
            ++dimRowCol;
            ++dimRowColM1;
            ++dimRowColP1;
            s1 = *(unsigned long*)(src + dimRowCol - 1);
            ss2 = *(unsigned long*)(src + dimRowCol);
            ss3 = *(unsigned long*)(src + dimRowCol + 1);
            s2 = ss2;
            s3 = ss3;
            srcR1 = (unsigned char)s2;
            s2 >>= 16;
            srcG1 = (unsigned char)s2;
            s2 >>= 16;
            srcB1 = (unsigned char)s2;
            r = -(unsigned char)s1 - srcR1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            g = -(unsigned char)s1 - srcG1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            b = -(unsigned char)s1 - srcB1 - (unsigned char)s3;

            dst[dimRowColM1].blue += b;
            dst[dimRowColM1].red += r;
            dst[dimRowColM1].green += g;

            //check if the value between 0 and 255
            dstB = dst[dimRowColM1].blue;
            dstR = dst[dimRowColM1].red;
            dstG = dst[dimRowColM1].green;

            if (dstB < 0)
                dst[dimRowColM1].blue = 0;
            else if (dstB > 255)
                dst[dimRowColM1].blue = 255;
            if (dstR < 0)
                dst[dimRowColM1].red = 0;
            else if (dstR > 255)
                dst[dimRowColM1].red = 255;
            if (dstG < 0)
                dst[dimRowColM1].green = 0;
            else if (dstG > 255)
                dst[dimRowColM1].green = 255;

            //add the value + 10* the same
            dst[dimRowCol].blue += 10*srcB1 + b;
            dst[dimRowCol].red += 10*srcR1 + r;
            dst[dimRowCol].green += 10*srcG1 + g;

            //initialize also
            d = ((((long)(unsigned short)(b)) << 16) << 16) +
            ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColP1)) = d;
            /////////////////////////////////////////////////////////////////////
            ++dimRowCol;
            ++dimRowColM1;
            ++dimRowColP1;
            ss1 = ss2;
            ss2 = ss3;
            ss3 = *(unsigned long*)(src + dimRowCol + 1);
            s1 = ss1;
            s2 = ss2;
            s3 = ss3;
            srcR1 = (unsigned char)s2;
            s2 >>= 16;
            srcG1 = (unsigned char)s2;
            s2 >>= 16;
            srcB1 = (unsigned char)s2;
            r = -(unsigned char)s1 - srcR1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            g = -(unsigned char)s1 - srcG1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            b = -(unsigned char)s1 - srcB1 - (unsigned char)s3;

            dst[dimRowColM1].blue += b;
            dst[dimRowColM1].red += r;
            dst[dimRowColM1].green += g;

            //check if the value between 0 and 255
            dstB = dst[dimRowColM1].blue;
            dstR = dst[dimRowColM1].red;
            dstG = dst[dimRowColM1].green;

            if (dstB < 0)
                dst[dimRowColM1].blue = 0;
            else if (dstB > 255)
                dst[dimRowColM1].blue = 255;
            if (dstR < 0)
                dst[dimRowColM1].red = 0;
            else if (dstR > 255)
                dst[dimRowColM1].red = 255;
            if (dstG < 0)
                dst[dimRowColM1].green = 0;
            else if (dstG > 255)
                dst[dimRowColM1].green = 255;

            //add the value + 10* the same
            dst[dimRowCol].blue += 10*srcB1 + b;
            dst[dimRowCol].red += 10*srcR1 + r;
            dst[dimRowCol].green += 10*srcG1 + g;

            //initialize also
            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColP1)) = d;
            ////////////////////////////////////////////////
            ++dimRowCol;
            ++dimRowColM1;
            ++dimRowColP1;
            ss1 = ss2;
            ss2 = ss3;
            ss3 = *(unsigned long*)(src + dimRowCol + 1);
            s1 = ss1;
            s2 = ss2;
            s3 = ss3;
            srcR1 = (unsigned char)s2;
            s2 >>= 16;
            srcG1 = (unsigned char)s2;
            s2 >>= 16;
            srcB1 = (unsigned char)s2;
            r = -(unsigned char)s1 - srcR1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            g = -(unsigned char)s1 - srcG1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            b = -(unsigned char)s1 - srcB1 - (unsigned char)s3;

            dst[dimRowColM1].blue += b;
            dst[dimRowColM1].red += r;
            dst[dimRowColM1].green += g;

            //check if the value between 0 and 255
            dstB = dst[dimRowColM1].blue;
            dstR = dst[dimRowColM1].red;
            dstG = dst[dimRowColM1].green;

            if (dstB < 0)
                dst[dimRowColM1].blue = 0;
            else if (dstB > 255)
                dst[dimRowColM1].blue = 255;
            if (dstR < 0)
                dst[dimRowColM1].red = 0;
            else if (dstR > 255)
                dst[dimRowColM1].red = 255;
            if (dstG < 0)
                dst[dimRowColM1].green = 0;
            else if (dstG > 255)
                dst[dimRowColM1].green = 255;

            //add the value + 10* the same
            dst[dimRowCol].blue += 10*srcB1 + b;
            dst[dimRowCol].red += 10*srcR1 + r;
            dst[dimRowCol].green += 10*srcG1 + g;

            //initialize also
            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColP1)) = d;
            ////////////////////////////////////////////
            ++dimRowCol;
            ++dimRowColM1;
            ++dimRowColP1;
            ss1 = ss2;
            ss2 = ss3;
            s3 = *(unsigned long*)(src + dimRowCol + 1);
            s1 = ss1;
            s2 = ss2;
            srcR1 = (unsigned char)s2;
            s2 >>= 16;
            srcG1 = (unsigned char)s2;
            s2 >>= 16;
            srcB1 = (unsigned char)s2;
            r = -(unsigned char)s1 - srcR1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            g = -(unsigned char)s1 - srcG1 - (unsigned char)s3;
            s1 >>= 16;
            s3 >>= 16;
            b = -(unsigned char)s1 - srcB1 - (unsigned char)s3;

            dst[dimRowColM1].blue += b;
            dst[dimRowColM1].red += r;
            dst[dimRowColM1].green += g;

            //check if the value between 0 and 255
            dstB = dst[dimRowColM1].blue;
            dstR = dst[dimRowColM1].red;
            dstG = dst[dimRowColM1].green;

            if (dstB < 0)
                dst[dimRowColM1].blue = 0;
            else if (dstB > 255)
                dst[dimRowColM1].blue = 255;
            if (dstR < 0)
                dst[dimRowColM1].red = 0;
            else if (dstR > 255)
                dst[dimRowColM1].red = 255;
            if (dstG < 0)
                dst[dimRowColM1].green = 0;
            else if (dstG > 255)
                dst[dimRowColM1].green = 255;

            //add the value + 10* the same
            dst[dimRowCol].blue += 10*srcB1 + b;
            dst[dimRowCol].red += 10*srcR1 + r;
            dst[dimRowCol].green += 10*srcG1 + g;

            //initialize also
            d = ((((long)(unsigned short)(b)) << 16) << 16) +
                ((long)(unsigned short)(g) << 16) + ((long)(unsigned short)(r));
            *((long *)(dst + dimRowColP1)) = d;
        }
        dimRowCol += 2 ;
        dimRowColM1 += 2;
        dimRowColP1 += 2;
    }
    //row = dim -2;
    dimRow = (dim - 2)*dim;
    dimRowM1 = dimRow - dim;
    for(col = 1; col < dim -1; ++col){
        b = src[dimRow + col -1].blue + src[dimRow + col].blue + src[dimRow + col + 1].blue;
        r = src[dimRow + col -1].red + src[dimRow + col].red + src[dimRow + col + 1].red;
        g = src[dimRow + col -1].green + src[dimRow + col].green + src[dimRow + col + 1].green;
        //add the value
        dst[dimRowM1 + col].blue -= b;
        dst[dimRowM1 + col].red -= r;
        dst[dimRowM1 + col].green -= g;

        //close
        dstB = dst[dimRowM1 + col].blue;
        dstR = dst[dimRowM1 + col].red;
        dstG = dst[dimRowM1 + col].green;

        if (dstB < 0)
            dst[dimRowM1 + col].blue = 0;
        else if (dstB > 255)
            dst[dimRowM1 + col].blue = 255;

        if (dstR < 0)
            dst[dimRowM1 + col].red = 0;
        else if (dstR > 255)
            dst[dimRowM1 + col].red = 255;

        if (dstG < 0)
            dst[dimRowM1 + col].green = 0;
        else if (dstG > 255)
            dst[dimRowM1 + col].green = 255;

        //add the value
        dst[dimRow + col].blue += -b + 10 * src[dimRow + col].blue;
        dst[dimRow + col].red += -r + 10 * src[dimRow + col].red;
        dst[dimRow + col].green += -g + 10 * src[dimRow + col].green;
    }

    //row = dim - 1
    dimRow +=dim;
    dimRowM1 +=dim;
    for(col = 1; col < dim -1; ++col){
        dst[dimRowM1 + col].blue -= src[dimRow + col - 1].blue + src[dimRow+ col].blue + src[dimRow + col + 1].blue;
        dst[dimRowM1 + col].red -= src[dimRow + col - 1].red + src[dimRow + col].red + src[dimRow + col + 1].red;
        dst[dimRowM1 + col].green -= src[dimRow + col - 1].green + src[dimRow + col].green + src[dimRow + col + 1].green;

        dstB = dst[dimRowM1 + col].blue;
        dstG = dst[dimRowM1 + col].green;
        dstR = dst[dimRowM1 + col].red;

        if (dstB < 0)
            dst[dimRowM1 + col].blue = 0;
        else if (dstB > 255)
            dst[dimRowM1 + col].blue = 255;

        if (dstR < 0)
            dst[dimRowM1 + col].red = 0;
        else if (dstR > 255)
            dst[dimRowM1 + col].red = 255;

        if (dstG < 0)
            dst[dimRowM1 + col].green = 0;
        else if (dstG > 255)
            dst[dimRowM1 + col].green = 255;
    }
}

void charsToPixels(pixel* pixels, int mn) {
    register unsigned long x;
    pixel* p = pixels;
    int mod = mn & 1;
    int index = 0, threeIndex = 0;
    char* imageData = image->data;
    if (mod = 1) {
        (*p).red = (unsigned char)imageData[threeIndex];
        (*p).green = (unsigned char)imageData[threeIndex + 1];
        (*p).blue = (unsigned char)imageData[threeIndex + 2];
        ++p;
        imageData += 3;
    }
    for (; index < mn ; index += 2, imageData += 6, ++p) {
        x = *((unsigned long*)(imageData));

        (*p).red = (unsigned char)x;
        x >>= 8;
        (*p).green = (unsigned char)x;
        x >>= 8;
        (*p).blue = (unsigned char)x;
        x >>= 8;
        (*(++p)).red = (unsigned char)x;
        x >>= 8;
        (*p).green = (unsigned char)x;
        x >>= 8;
        (*p).blue = (unsigned char)x;
    }
}

//like chars to pixel but also calculate the sim value of each pixel
void charsToPixels7(Image *charsImg, pixel* pixels, int mn) {
    register unsigned long x;
    pixel* p = pixels;
    int mod = mn & 1;
    int index = 0, threeIndex = 0;
    char* imageData = image->data;
    if (mod = 1) {
        (*p).red = (unsigned char)imageData[threeIndex];
        (*p).green = (unsigned char)imageData[threeIndex + 1];
        (*p).blue = (unsigned char)imageData[threeIndex + 2];
        (*p).sum = (*p).blue + (*p).green + (*p).red;
        ++p;
        imageData += 3;
    }
    for (; index < mn ; index += 2, imageData += 6, ++p) {
        x = *((unsigned long*)(imageData));

        (*p).red = (unsigned char)x;
        x >>= 8;
        (*p).green = (unsigned char)x;
        x >>= 8;
        (*p).blue = (unsigned char)x;
        (*p).sum = (*p).blue + (*p).green + (*p).red;
        x >>= 8;
        (*(++p)).red = (unsigned char)x;
        x >>= 8;
        (*p).green = (unsigned char)x;
        x >>= 8;
        (*p).blue = (unsigned char)x;
        (*p).sum = (*p).blue + (*p).green + (*p).red;
    }
}

void pixelsToChars(pixel* pixels, int mn) {
    register unsigned int x1;
    register unsigned int x2;
    register unsigned int x3;
    register unsigned int x4;
    pixel * p = pixels;
    int mod = mn % 4;
    int index = 0;
    char* imageData = image->data;
    imageData -= 3;
    --p;
    for(;index < mod; ++index) {
        ++p;
        x1 = ((*(p)).blue << 16) + ((*(p)).green << 8) + ((*(p)).red);
        imageData += 3;
        *((int *)imageData) = x1;
    }
    for (; index < mn ; index += 4) {
        //lop1
        ++p;
        x1 = ((*(p)).blue << 16) + ((*(p)).green << 8) + ((*(p)).red);
        imageData += 3;
        *((int *)imageData) = x1;

        ++p;
        x2 = ((*(p)).blue << 16) + ((*(p)).green << 8) + ((*(p)).red);
        imageData += 3;
        *((int *)imageData) = x2;

        ++p;
        x3 = ((*(p)).blue << 16) + ((*(p)).green << 8) + ((*(p)).red);
        imageData += 3;
        *((int *)imageData) = x3;

        ++p;
        x4 = ((*(p)).blue << 16) + ((*(p)).green << 8) + ((*(p)).red);
        imageData += 3;
        *((int *)imageData) = x4;
    }
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {
    int mn = m*n;
    int size = mn*sizeof(pixel);
    pixel* pixelsImg = (pixel*)malloc(2*size);
    pixel* backupOrg = pixelsImg + mn;

	if (flag == '1') {

        charsToPixels(backupOrg, mn);

        blurKernelSmoothWithoutFilter9( m, backupOrg, pixelsImg);
        sharpKernelSmoothWithoutFilter( m ,pixelsImg, backupOrg);

        pixelsToChars(pixelsImg, mn);
		writeBMP(image, srcImgpName, blurRsltImgName);
        pixelsToChars(backupOrg, mn);
		writeBMP(image, srcImgpName, sharpRsltImgName);

	} else {

        charsToPixels7(image, backupOrg, mn);

        blurKernelSmoothWithFilter77( m, backupOrg, pixelsImg);
        sharpKernelSmoothWithoutFilter( m ,pixelsImg, backupOrg);

        pixelsToChars(pixelsImg, mn);
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);
        pixelsToChars(backupOrg, mn);
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);	
	}

    free(pixelsImg);
}

