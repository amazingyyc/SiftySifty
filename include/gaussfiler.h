/**
 * Created by yanyuanchi on 2017/3/21.
 */
#ifndef SIFTYSIFTY_GAUSSFILER_H
#define SIFTYSIFTY_GAUSSFILER_H

namespace SiftySifty {

/**
 * use 3 box filter to fitting gauss filter
 */
void gaussFilterBy3BoxFilter(short *src, short *dst, int width, int height, float sigma);

/**
 * usr IIR filter to fitting gauss filter
 */
void gaussFilterByIIRFilter(short *src, short *dst, int width, int height, float sigma);

}

#endif //SIFTYSIFTY_GAUSSFILER_H
