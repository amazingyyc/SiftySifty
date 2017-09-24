/**
 * Created by yanyuanchi on 2017/3/23.
 */

#ifndef SIFTYSIFTY_IIRFILTER_H
#define SIFTYSIFTY_IIRFILTER_H

namespace SiftySifty {

void IIRFilter(int16_t *src, int16_t *dst, int width, int height, float sigma);

void IIRFilter(uint8_t *src, uint8_t *dst, int width, int height, float sigma);

}

#endif //SIFTYSIFTY_IIRFILTER_H
