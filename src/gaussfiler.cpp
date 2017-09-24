/**
 * Created by yanyuanchi on 2017/3/21.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "boxfilter.h"
#include "iirfilter.h"
#include "gaussfiler.h"

namespace SiftySifty {

/**
 * use 3 box filter to fitting gauss filter
 */
void gaussFilterBy3BoxFilter(int16_t *src, int16_t *dst, int width, int height, float sigma) {
    /**
     * get the radius for 3 box filter
     * ref:http://blog.ivank.net/fastest-gaussian-blur.html
     */
    float wIdeal = sqrt(12.0 * sigma * sigma / 3 + 1.0);
    int wl = floor(wIdeal);
    
    if (0 == wl % 2) {
        wl--;
    }
    
    int wu = wl + 2;
    
    float mIdeal = (12.0 * sigma * sigma - 3 * wl * wl - 4 * 3 * wl - 3 * 3) / (-4 * wl - 4);
    int m = round(mIdeal);
    
    int radius[3];
    for (int i = 0; i < 3; ++i) {
        radius[i] = (i < m ? wl : wu) / 2;
    }
    
    short *tmp = (short *) malloc(sizeof(short) * width * height);
    
    boxFilter(src, dst, width, height, radius[0]);
    boxFilter(dst, tmp, width, height, radius[1]);
    boxFilter(tmp, dst, width, height, radius[2]);
    
    free(tmp);
}

/**
 * use the IIR method to fit the gauss filter
 * @param src
 * @param dst
 * @param width
 * @param height
 * @param sigma
 */
void gaussFilterByIIRFilter(short *src, short *dst, int width, int height, float sigma) {
    IIRFilter(src, dst, width, height, sigma);
}

}