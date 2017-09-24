/**
 * Created by yanyuanchi on 2017/1/18.
 */
#ifndef SIFT_LINEARFILTER_H
#define SIFT_LINEARFILTER_H

namespace SiftySifty {
/**
 * filter the rows
 */
void linearFilterHorizon(uint8_t *src,
                         uint8_t *dst,
                         int width,
                         int height,
                         int (*mult)[256],
                         int delta,
                         int shift,
                         int size);

/**
 * filter the cols
 */
void linearFilterVertical(uint8_t *src,
                          uint8_t *dst,
                          int width,
                          int height,
                          int (*mult)[256],
                          int delta,
                          int shift,
                          int size);
}
#endif //SIFT_LINEARFILTER_H
