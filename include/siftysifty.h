/**
 * Created by yanyuanchi on 2017/3/22.
 */
#ifndef SIFTYSIFTY_SIFTYSIFTY_H
#define SIFTYSIFTY_SIFTYSIFTY_H

#include <iostream>
#include <vector>

#include "utils.h"
#include "structs.h"
#include "gaussfiler.h"
#include "siftysifty.h"

using namespace std;

namespace SiftySifty {
/**the gray image will be scale by the SIFT_IMAGE_SCALE*/
static const int SIFT_IMAGE_SCALE_SHIFT = 6;
static const int SIFT_IMAGE_SCALE = (1 << SIFT_IMAGE_SCALE_SHIFT);

/**默认的金字塔每一层需要计算特征点的图片个数*/
/**the pic's number on one layer*/
static const int SIFT_OCTAVE_LAYERS = 3;

/**the sigma of sift*/
static const float SIFT_SIGMA = 1.6f;

/**the base pic's sigma*/
static const float SIFT_INIT_SIGMA = 0.5f;

/**sift's contrast threshold*/
static const float SIFT_CONTRAST_THRESHOLD = 0.04f;

/**the edge threshold*/
static const float SIFT_EDGE_THESHOLD = 10;

/**if the inited pic will be doubled*/
static const bool SIFT_DOUBLE_INITED_IMAGE = true;

/**the region width of for descriptor*/
static const int SIFT_DESCRIPTOR_WIDTH = 4;

/**the number of image region*/
static const int SIFT_DESCRIPTOR_HIST_BIN = 8;

static const float SIFT_ORIENTATION_PEAK_RATIO = 0.8f;

/**360 splited to 36*/
static const int SIFT_ORIENTATION_HIST_BINS = 36;

/**adjust 5 time*/
static const int SIFT_MAX_ADJUST_STEP = 5;

/**the border of the image*/
static const int SIFT_IMAGE_BORDER = 5;

static const float SIFT_ORIENTATION_SIGMA_FCTER = 1.5f;

/**the radius of hist SIFT_ORIENTATION_RADIUS * sigma*/
static const float SIFT_ORIENTATION_RADIUS = (3 * SIFT_ORIENTATION_SIGMA_FCTER);

/**SIFT_DESCRPTOR_SCAE_FCTER * sigma*/
static const float SIFT_DESCRIPTOR_SCAE_FCTER = 3.0f;

static const float SIFT_DESCRIPTOR_MAGNITUDE_THRESHOLD = 0.2f;

static const float SIFT_DESCRIPTOR_FCTOR = 512.0f;

void sift(uint8_t *image,
          int width, int height,
          vector<SiftySifty::KeyPoint> &keyPoints,
          int octaveLayers,
          float sigma,
          float contrastThreshold,
          int edgeThreshold,
          bool doubleInitImage,
          int descriptorWidth,
          int descriptorHistBin);

void sift(uint8_t *image, int width, int height, vector<SiftySifty::KeyPoint> &keyPoints);

void initSpeed(uint8_t *src, int width, int height);

}

#endif //SIFTYSIFTY_SIFTYSIFTY_H



















