/**
 * Created by yanyuanchi on 2017/3/21.
 */
#ifndef SIFTYSIFTY_BOXFILTER_H
#define SIFTYSIFTY_BOXFILTER_H

namespace SiftySifty {
/**
 * box blur
 */
void boxFilter(int16_t *src, int16_t *dst, int width, int height, int radius);

}

#endif //SIFTYSIFTY_BOXFILTER_H
