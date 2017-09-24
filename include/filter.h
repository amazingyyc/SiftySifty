/**
 * Created by yanyuanchi on 2017/3/23.
 */
#ifndef SIFTYSIFTY_FILTER_H
#define SIFTYSIFTY_FILTER_H

namespace SiftySifty {

static const int FILTER_SHIFT = 16;
static const int FILTER_SCALE = (1 << FILTER_SHIFT);
static const int FILTER_DELTA = (1 << (FILTER_SHIFT - 1));

}

#endif //SIFTYSIFTY_FILTER_H
