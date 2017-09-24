/**
 * Created by yanyuanchi on 2017/3/22.
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <climits>
#include <algorithm>

#include "utils.h"
#include "gaussfiler.h"
#include "imageutils.h"
#include "iirfilter.h"
#include "siftysifty.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace SiftySifty {

/**
 * gauss to Mat<short>
 */
void gaussFilter(Mat<int16_t> *src, Mat<short> *dst, float sigma) {
    if (nullptr == src || nullptr == dst || 0 >= sigma || src->width != dst->width || src->height != dst->height) {
        return;
    }
    
    gaussFilterByIIRFilter(src->data, dst->data, src->width, src->height, sigma);
}

/**
 * gauss blur on the mat
 */
void gaussFilter(Mat<int16_t> *mat, float sigma) {
    Mat<int16_t> *tmp = newMat<int16_t>(mat->width, mat->height);
    
    memcpy(tmp->data, mat->data, sizeof(int16_t) * mat->width * mat->height);
    
    gaussFilter(tmp, mat, sigma);
    
    deleteMat<int16_t>(tmp);
}

/**
 * create the base image of sift
 * @param src the gray image with uint8_t
 * @param width width
 * @param height height
 * @param doubleImage if double the image
 * @param sigma the init sigma
 * @param shift the base image will be src * (1 << shift)
 * @return
 */
Mat<int16_t> *initBaseImage(uint8_t *src, int width, int height, bool doubleImage, double sigma, int shift) {
    Mat<int16_t> *base = newMat<int16_t>(width, height);
    
    /**scale the src*/
    scaleMatByShift<uint8_t, int16_t>(src, base->data, width, height, shift);
    
    if (doubleImage) {
        float diffSigma = (float) sqrt(sigma * sigma - 4.0 * SIFT_INIT_SIGMA * SIFT_INIT_SIGMA);
        
        Mat<int16_t> *doubleBase = newMat<int16_t>(base->width * 2, base->height * 2);
        
        resizeMat2<int16_t>(base, doubleBase);
        gaussFilter(doubleBase, diffSigma);
        
        deleteMat<int16_t>(base);
        
        return doubleBase;
    } else {
        float diffSigma = (float) sqrt(sigma * sigma - SIFT_INIT_SIGMA * SIFT_SIGMA);
        
        gaussFilter(base, diffSigma);
        
        return base;
    }
}

vector<vector<SiftySifty::Mat<int16_t> *> > buildGaussPyramid(SiftySifty::Mat<int16_t> *base,
                                                              const int octave,
                                                              const int octaveLayers,
                                                              const float sigma) {
    double sig[octaveLayers + 3];
    sig[0] = sigma;
    
    double k = pow(2.0, 1.0 / octaveLayers);
    for (int i = 1; i < octaveLayers + 3; i++) {
        double sigPrev = pow(k, (double) (i - 1)) * sigma;
        double sigTotal = sigPrev * k;
        
        sig[i] = sqrt(sigTotal * sigTotal - sigPrev * sigPrev);
    }
    
    vector<vector<SiftySifty::Mat<int16_t> *> > gaussPyramid(octave,
                                                             vector<SiftySifty::Mat<int16_t> *>(octaveLayers + 3));
    
    for (int y = 0; y < octave; ++y) {
        for (int x = 0; x < (octaveLayers + 3); ++x) {
            if (0 == y && 0 == x) {
                gaussPyramid[y][x] = base;
            } else if (0 == x) {
                Mat<int16_t> *pre = gaussPyramid[y - 1][octaveLayers];
                Mat<int16_t> *mat = newMat<int16_t>(pre->width >> 1, pre->height >> 1);
                
                halfSampleMat<int16_t>(pre, mat);
                
                gaussPyramid[y][x] = mat;
            } else {
                Mat<int16_t> *pre = gaussPyramid[y][x - 1];
                Mat<int16_t> *mat = newMat<int16_t>(pre->width, pre->height);
                
                gaussFilter(pre, mat, sig[x]);
                
                gaussPyramid[y][x] = mat;
            }
        }
    }
    
    return gaussPyramid;
}

/**
 * build dogPyramid
 */
vector<vector<SiftySifty::Mat<int16_t> *> > buildDoGPyramid(vector<vector<SiftySifty::Mat<int16_t> *> > &gaussPyramid,
                                                            const int octave,
                                                            const int octaveLayers) {
    vector<vector<SiftySifty::Mat<int16_t> *> > doGPyramid(octave,
                                                           vector<SiftySifty::Mat<int16_t> *>(octaveLayers + 2));
    
    for (int y = 0; y < octave; ++y) {
        for (int x = 0; x < (octaveLayers + 2); ++x) {
            SiftySifty::Mat<int16_t> *src1 = gaussPyramid[y][x + 1];
            SiftySifty::Mat<int16_t> *src2 = gaussPyramid[y][x];
            
            SiftySifty::Mat<int16_t> *dst = newMat<int16_t>(src1->width, src1->height);
            
            subMat<int16_t>(src1->data, src2->data, dst->data, src1->width, src1->height);
            
            doGPyramid[y][x] = dst;
        }
    }
    
    return doGPyramid;
}

/**
 * adjust the real extrema point
 * @param doGPyramid dog pyramid
 * @param keyPoint the current keypoint
 * @param r the row of keypoint
 * @param c the col of keypoint
 * @param layer the layer of the keypoint
 * @param octaveLayers the total images of the layer
 * @param curOctave the current octave
 * @param contrastThreshold contrast threshold
 * @param edgeThreshold edge threshold
 * @param sigma
 * @param offset
 * @return true: is a keypint false:not
 */
bool adjustExtremaPoint(vector<vector<SiftySifty::Mat<int16_t> * > > &doGPyramid,
                        KeyPoint &keyPoint,
                        int &r, int &c, int &l,
                        int octaveLayers,
                        int curOctave,
                        float contrastThreshold,
                        float edgeThreshold,
                        float sigma,
                        int *offset) {
    /**
     * the origin sift paper use the float to store the image,
     * so the imageScale have to be scaled by the SIFT_IMAGE_SCALE
     */
    const float imageScale = 1.0f / 255.0f * SIFT_IMAGE_SCALE;
    
    /**benn used to calcualte the first-order derivative*/
    const float deriveScale = 0.5f * imageScale;
    
    /**second-order derivative*/
    const float secondDeriveScale = imageScale;
    
    /**cross--order derivative*/
    const float crossDeriveScale = 0.25f * imageScale;
    
    float xR = 0, xC = 0, xL = 0;
    
    bool confirm = false;
    
    for (int i = 0; i < SIFT_MAX_ADJUST_STEP; ++i) {
        int width = doGPyramid[curOctave][l]->width;
        int height = doGPyramid[curOctave][l]->height;
        
        int16_t *cur = doGPyramid[curOctave][l]->data + width * r + c;
        int16_t *pre = doGPyramid[curOctave][l - 1]->data + width * r + c;
        int16_t *nex = doGPyramid[curOctave][l + 1]->data + width * r + c;
        
        /**
         * calculate the derive of x, y, sigma
         *         |
         *         |
         *         |
         * ---------------->x
         *         |
         *         |
         *         |
         *         ^ y
         */
        float dx = (cur[offset[5]] - cur[offset[3]]) * deriveScale;
        float dy = (cur[offset[7]] - cur[offset[1]]) * deriveScale;
        float ds = (nex[0] - pre[0]) * deriveScale;
        
        /**
         * calculate dxx dxy dxs dyx dyy dys dsx dsy dss
         */
        float value2 = 2.0f * cur[0];
        float dxx = (cur[offset[5]] + cur[offset[3]] - value2) * secondDeriveScale;
        float dyy = (cur[offset[7]] + cur[offset[1]] - value2) * secondDeriveScale;
        float dss = (nex[0] + pre[0] - value2) * secondDeriveScale;
        
        float dxy = (cur[offset[8]] + cur[offset[0]] - cur[offset[2]] - cur[offset[6]]) * crossDeriveScale;
        float dxs = (nex[offset[5]] + pre[offset[3]] - nex[offset[3]] - pre[offset[5]]) * crossDeriveScale;
        float dys = (nex[offset[7]] + pre[offset[1]] - nex[offset[1]] - pre[offset[7]]) * crossDeriveScale;
        
        /**
         * X = -[dx, dy, ds] ^ T * ([Dxx, Dxy, Dxs]) ^ -1
         *                         ([Dyx, Dyy, Dys])
         *                         ([Dsx, Dsy, Dss])
         */
        float detD = dxx * dyy * dss + dxy * dys * dxs + dxs * dxy * dys -
                     dxx * dys * dys - dxy * dxy * dss - dxs * dyy * dxs;
        
        if (fabsf(detD) < 1e-6) {
            return false;
        }
        
        detD = 1.0 / detD;
        
        xC = dx * (dyy * dss - dys * dys) + dy * (-dxy * dss + dys * dxs) + ds * (dxy * dys - dxs * dyy);
        xR = dx * (dys * dxs - dss * dxy) + dy * (-dxs * dxs + dss * dxx) + ds * (dxs * dxy - dxx * dys);
        xL = dx * (dxy * dys - dxs * dyy) + dy * (-dxx * dys + dxs * dxy) + ds * (dxx * dyy - dxy * dxy);
        
        xC = -detD * xC;
        xR = -detD * xR;
        xL = -detD * xL;
        
        if (fabsf(xC) < 0.5f && fabsf(xR) < 0.5f && fabsf(xL) < 0.5f) {
            confirm = true;
            break;
        }
        
        if (std::abs(xC) > (float) (INT_MAX / 3) ||
            std::abs(xR) > (float) (INT_MAX / 3) ||
            std::abs(xL) > (float) (INT_MAX / 3)) {
            return false;
        }
        
        r += (int) (roundf(xR));
        c += (int) (roundf(xC));
        l += (int) (roundf(xL));
        
        /**out of border*/
        if (l < 1 || l > octaveLayers
            || c < SIFT_IMAGE_BORDER || c >= width - SIFT_IMAGE_BORDER
            || r < SIFT_IMAGE_BORDER || r >= height - SIFT_IMAGE_BORDER) {
            return false;
        }
    }
    
    if (!confirm) {
        return false;
    }
    
    int width = doGPyramid[curOctave][l]->width;
    
    int16_t *cur = doGPyramid[curOctave][l]->data + width * r + c;
    int16_t *pre = doGPyramid[curOctave][l - 1]->data + width * r + c;
    int16_t *nex = doGPyramid[curOctave][l + 1]->data + width * r + c;
    
    float dx = (cur[offset[5]] - cur[offset[3]]) * deriveScale;
    float dy = (cur[offset[7]] - cur[offset[1]]) * deriveScale;
    float ds = (nex[0] - pre[0]) * deriveScale;
    
    float response = cur[0] * imageScale + 0.5f * (dx * xC + dy * xR + ds * xL);
    
    if (fabsf(response) * octaveLayers < contrastThreshold) {
        return false;
    }
    
    float value2 = 2.0f * cur[0];
    float dxx = (cur[offset[5]] + cur[offset[3]] - value2) * secondDeriveScale;
    float dyy = (cur[offset[7]] + cur[offset[1]] - value2) * secondDeriveScale;
    float dxy = (cur[offset[8]] + cur[offset[0]] - cur[offset[2]] - cur[offset[6]]) * crossDeriveScale;
    
    float tr = dxx + dyy;
    float det = dxx * dyy - dxy * dxy;
    
    if (0 > det || tr * tr * edgeThreshold >= (edgeThreshold + 1) * (edgeThreshold + 1) * det) {
        return false;
    }
    
    keyPoint.x = (c + xC) * (1 << curOctave);
    keyPoint.y = (r + xR) * (1 << curOctave);
    
    keyPoint.octaveX = c;
    keyPoint.octaveY = r;
    
    keyPoint.octave = curOctave;
    keyPoint.octaveLayer = l;
    
    keyPoint.octaveLayersShift = xL;
    
    keyPoint.size = sigma * powf(2.f, (l + xL) / octaveLayers) * (1 << curOctave);
    keyPoint.octaveSize = sigma * powf(2.f, (l + xL) / octaveLayers);
    
    /**[-1, 1]*/
    keyPoint.response = response;
    
    return true;
}

/**
 * calculate hist
 */
float calculateHist(SiftySifty::Mat<int16_t> *image, int x, int y, int radius, float sigma, float *hist, int n) {
    int16_t *imageData = image->data;
    int width = image->width;
    int height = image->height;
    
    int length = (2 * radius + 1) * (2 * radius + 1);
    
    float scale = -1.0f / (2.0f * sigma * sigma);
    
    float *buffer = (float *) malloc(sizeof(float) * (5 * length + n + 4));
    
    float *weights = buffer, *dx = weights + length, *dy = dx + length, *ori = dy + length, *mag = ori + length;
    float *tmpHist = mag + length + 2;
    
    memset(tmpHist, 0, sizeof(float) * n);
    
    int r, c;
    int realLength = 0;
    for (int i = -radius; i <= radius; ++i) {
        r = y + i;
        if (0 >= r || r >= (height - 1)) {
            continue;
        }
        
        for (int j = -radius; j <= radius; ++j) {
            c = x + j;
            if (0 >= c || c >= (width - 1)) {
                continue;
            }
            /**
             *             ^ y
             *             |
             *             |
             *             |
             * -----------------------> x
             *             |
             *             |
             *             |
             */
            dx[realLength] = (float) (imageData[r * width + c + 1] - imageData[r * width + c - 1]);
            dy[realLength] = (float) (imageData[(r - 1) * width + c] - imageData[(r + 1) * width + c]);
            
            weights[realLength] = (i * i + j * j) * scale;
            
            realLength++;
        }
    }
    
    for (int i = 0; i < realLength; ++i) {
        weights[i] = expf(weights[i]);
        mag[i] = sqrtf(dx[i] * dx[i] + dy[i] * dy[i]);
        ori[i] = atan2f360(dy[i], dx[i]);
    }
    
    for (int i = 0; i < realLength; ++i) {
        int index = static_cast<int>(roundf(ori[i] * n / 360.0f));
        
        if (index > n) {
            index -= n;
        }
        
        if (0 > index) {
            index += n;
        }
        
        tmpHist[index] += (weights[i] * mag[i]);
    }
    
    tmpHist[-1] = tmpHist[n - 1];
    tmpHist[-2] = tmpHist[n - 2];
    tmpHist[n] = tmpHist[0];
    tmpHist[n + 1] = tmpHist[1];
    
    for (int i = 0; i < n; ++i) {
        hist[i] = (tmpHist[i - 2] + tmpHist[i + 2]) * (1.f / 16.f) +
                  (tmpHist[i - 1] + tmpHist[i + 1]) * (4.f / 16.f) + tmpHist[i] * (6.f / 16.f);
    }
    
    float maxValue = hist[0];
    for (int i = 0; i < n; ++i) {
        maxValue = max_value(hist[i], maxValue);
    }
    
    free(buffer);
    
    return maxValue;
}

void findKeyPoints(vector<KeyPoint> &kpts, KeyPoint &keyPoint, float threshold, float *hist, int n) {
    for (int i = 0; i < n; ++i) {
        int left = (i > 0) ? (i - 1) : (n - 1);
        int right = (i < (n - 1)) ? (i + 1) : 0;
        
        if (hist[i] > hist[left] && hist[i] >= hist[right] && hist[i] >= threshold) {
            float bin = i + 0.5f * (hist[left] - hist[right]) / (hist[left] - 2 * hist[i] + hist[right]);
            
            bin = (bin < 0) ? (n + bin) : ((bin >= n) ? (bin - n) : bin);
            
            keyPoint.angle = 360.0f - (360.0f / n) * bin;
            
            if (fabsf(keyPoint.angle - 360.0f) < FLT_EPSILON) {
                keyPoint.angle = 0.f;
            }
            
            kpts.push_back(keyPoint);
        }
    }
}

/**
 * find the extrema point on one pic
 * @param gaussPyramid gauss pyramid
 * @param doGPyramid dog pyramid
 * @param keyPoints store the keyPoints
 * @param octave the size of pyramid
 * @param octaveLayers the image's number in one layer
 * @param contrastThreshold the contrast threshold int sift paper
 * @param edgeThreshold the edge threshold in sift paper
 * @param valueThreshold the image value threshold
 * @param curOctave the current octave
 * @param curLayer the curent layer
 * @param n
 * @param sigma
 * @param offset
 */
void findExtremaPointOne(vector<vector<SiftySifty::Mat<int16_t> * > > &gaussPyramid,
                         vector<vector<SiftySifty::Mat<int16_t> * > > &doGPyramid,
                         vector<SiftySifty::KeyPoint> &keyPoints,
                         int octave,
                         int octaveLayers,
                         float contrastThreshold,
                         float edgeThreshold,
                         int valueThreshold,
                         int curOctave,
                         int curLayer,
                         int n,
                         float sigma,
                         int *offset) {
    /**get the pre current next image*/
    Mat<int16_t> *pre = doGPyramid[curOctave][curLayer - 1];
    Mat<int16_t> *cur = doGPyramid[curOctave][curLayer];
    Mat<int16_t> *nex = doGPyramid[curOctave][curLayer + 1];
    
    int width = cur->width;
    int height = cur->height;
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    /**remove the border*/
    int stride = max((int) (roundf(1.0f * (height - SIFT_IMAGE_BORDER - SIFT_IMAGE_BORDER) / maxThreadNum)), 1);

#ifdef _OPENMP
    omp_lock_t lock;
    omp_init_lock(&lock);
#endif

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride + SIFT_IMAGE_BORDER;
        int end = (threadIndex == (maxThreadNum - 1)) ?
                  (height - SIFT_IMAGE_BORDER) : min_value(start + stride, height - SIFT_IMAGE_BORDER);
        
        float hist[n];
        
        KeyPoint keyPoint;
        vector<KeyPoint> kpts;
        
        int16_t *prePtr = pre->data + start * width + SIFT_IMAGE_BORDER;
        int16_t *curPtr = cur->data + start * width + SIFT_IMAGE_BORDER;
        int16_t *nexPtr = nex->data + start * width + SIFT_IMAGE_BORDER;
        
        for (int y = start; y < end; ++y) {
            int16_t *preData = prePtr - 1;
            int16_t *curData = curPtr - 1;
            int16_t *nexData = nexPtr - 1;
            
            for (int x = SIFT_IMAGE_BORDER; x < (width - SIFT_IMAGE_BORDER); ++x) {
                preData++;
                curData++;
                nexData++;
                
                int val = curData[0];
                
                if (!(abs(val) > valueThreshold &&
                      ((val > 0 && val >= curData[offset[0]] && val >= curData[offset[1]] &&
                        val >= curData[offset[2]] && val >= curData[offset[3]] && val >= curData[offset[5]] &&
                        val >= curData[offset[6]] && val >= curData[offset[7]] && val >= curData[offset[8]] &&
                        val >= preData[offset[0]] && val >= preData[offset[1]] && val >= preData[offset[2]] &&
                        val >= preData[offset[3]] && val >= preData[offset[4]] && val >= preData[offset[5]] &&
                        val >= preData[offset[6]] && val >= preData[offset[7]] && val >= preData[offset[8]] &&
                        val >= nexData[offset[0]] && val >= nexData[offset[1]] && val >= nexData[offset[2]] &&
                        val >= nexData[offset[3]] && val >= nexData[offset[4]] && val >= nexData[offset[5]] &&
                        val >= nexData[offset[6]] && val >= nexData[offset[7]] && val >= nexData[offset[8]]) ||
                       (val < 0 && val <= curData[offset[0]] && val <= curData[offset[1]] &&
                        val <= curData[offset[2]] && val <= curData[offset[3]] && val <= curData[offset[5]] &&
                        val <= curData[offset[6]] && val <= curData[offset[7]] && val <= curData[offset[8]] &&
                        val <= preData[offset[0]] && val <= preData[offset[1]] && val <= preData[offset[2]] &&
                        val <= preData[offset[3]] && val <= preData[offset[4]] && val <= preData[offset[5]] &&
                        val <= preData[offset[6]] && val <= preData[offset[7]] && val <= preData[offset[8]] &&
                        val <= nexData[offset[0]] && val <= nexData[offset[1]] && val <= nexData[offset[2]] &&
                        val <= nexData[offset[3]] && val <= nexData[offset[4]] && val <= nexData[offset[5]] &&
                        val <= nexData[offset[6]] && val <= nexData[offset[7]] && val <= nexData[offset[8]])))) {
                    continue;
                }
                
                int extremaR = y, extremaC = x, extremaL = curLayer;
                
                if (!adjustExtremaPoint(doGPyramid,
                                        keyPoint,
                                        extremaR, extremaC, extremaL,
                                        octaveLayers,
                                        curOctave,
                                        contrastThreshold,
                                        edgeThreshold,
                                        sigma,
                                        offset)) {
                    continue;
                }
                
                float octaveSigma = keyPoint.octaveSize;
                
                float maxValue = calculateHist(gaussPyramid[curOctave][extremaL],
                                               extremaC,
                                               extremaR,
                                               (int) roundf(SIFT_ORIENTATION_SIGMA_FCTER * octaveSigma),
                                               SIFT_ORIENTATION_RADIUS * octaveSigma,
                                               hist,
                                               n);
                
                findKeyPoints(kpts, keyPoint, maxValue * SIFT_ORIENTATION_PEAK_RATIO, hist, n);
            }
            
            prePtr += width;
            curPtr += width;
            nexPtr += width;
        }
        
        if (!kpts.empty()) {
#ifdef _OPENMP
            omp_set_lock(&lock);
#endif
            for (auto iter = kpts.begin(); iter < kpts.end(); ++iter) {
                keyPoints.push_back(*iter);
            }

#ifdef _OPENMP
            omp_unset_lock(&lock);
#endif
        }
    }

#ifdef _OPENMP
    omp_destroy_lock(&lock);
#endif
}

/**
 * find the extrema point
 */
void findExtremaPoint(vector<vector<SiftySifty::Mat<int16_t> * >> &gaussPyramid,
                      vector<vector<SiftySifty::Mat<int16_t> * >> &doGPyramid,
                      vector<SiftySifty::KeyPoint> &keyPoints,
                      int octave,
                      int octaveLayers,
                      float contrastThreshold,
                      float edgeThreshold,
                      float sigma) {
    int n = SIFT_ORIENTATION_HIST_BINS;
    int threshold = (int) (0.5f * contrastThreshold / octaveLayers * 255 * SIFT_IMAGE_SCALE);
    
    int offset[9];
    
    for (int o = 0; o < octave; ++o) {
        int width = doGPyramid[o][0]->width;
        
        offset[0] = -width - 1;
        offset[1] = -width;
        offset[2] = -width + 1;
        offset[3] = -1;
        offset[4] = 0;
        offset[5] = 1;
        offset[6] = width - 1;
        offset[7] = width;
        offset[8] = width + 1;
        
        for (int l = 1; l <= octaveLayers; ++l) {
            findExtremaPointOne(gaussPyramid,
                                doGPyramid,
                                keyPoints,
                                octave,
                                octaveLayers,
                                contrastThreshold,
                                edgeThreshold,
                                threshold,
                                o,
                                l,
                                n,
                                sigma,
                                offset);
        }
    }
}

/**
 * if the inited image is been doubled, than resize the size
 */
void resizeKeyPoints(vector<SiftySifty::KeyPoint> &keyPoints) {
    auto iter = keyPoints.begin();
    
    for (; iter < keyPoints.end(); ++iter) {
        (*iter).x /= 2.0f;
        (*iter).y /= 2.0f;
        (*iter).size /= 2.0f;
    }
}

/**
 * sort the keypoint
 */
struct KeyPointCMP {
    vector<SiftySifty::KeyPoint> keyPoints;
    
    KeyPointCMP(const vector<SiftySifty::KeyPoint> &keyPoints) {
        this->keyPoints = keyPoints;
    }
    
    bool operator()(int i, int j) const {
        SiftySifty::KeyPoint kp1 = keyPoints[i];
        SiftySifty::KeyPoint kp2 = keyPoints[j];
        
        if (kp1.x != kp2.x) {
            return kp1.x < kp2.x;
        }
        
        if (kp1.y != kp2.y) {
            return kp1.y < kp2.y;
        }
        
        if (kp1.size != kp2.size) {
            return kp1.size < kp2.size;
        }
        
        if (kp1.angle != kp2.angle) {
            return kp1.angle < kp2.angle;
        }
        
        return i < j;
    }
};

/**
 * remove the doubled keypoint
 * @param keyPoints
 */
void removeDoubleKeyPoints(vector<SiftySifty::KeyPoint> &keyPoints) {
    int i, j;
    int n = keyPoints.size();
    
    /**sorted keypoint*/
    vector<int> sortedIndex(n);
    
    /**mark the keypoint*/
    vector<uint8_t> map(n, 1);
    
    for (int i = 0; i < sortedIndex.size(); ++i) {
        sortedIndex[i] = i;
    }
    
    /**sort the keypoint*/
    std::sort(sortedIndex.begin(), sortedIndex.end(), KeyPointCMP(keyPoints));
    
    for (i = 1, j = 0; i < n; ++i) {
        KeyPoint kp1 = keyPoints[sortedIndex[j]];
        KeyPoint kp2 = keyPoints[sortedIndex[i]];
        
        if (kp1.x != kp2.x
            || kp1.y != kp2.y
            || kp1.size != kp2.size
            || kp1.angle != kp2.angle) {
            j = i;
        } else {
            map[i] = 0;
        }
    }
    
    for (i = 0, j = 0; i < n; ++i) {
        if (1 == map[i]) {
            if (i != j) {
                keyPoints[j] = keyPoints[i];
            }
            
            j++;
        }
    }
    
    keyPoints.resize(j);
}

/**
 * calculate the keypoint's descriptor
 * @param image
 * @param x
 * @param y
 * @param angle
 * @param scale
 * @param d
 * @param n
 * @param descriptor
 */
void calculateDescriptorOne(Mat<int16_t> *image,
                            int x, int y,
                            float angle,
                            float scale,
                            int d,
                            int n,
                            float *descriptor) {
    int16_t *data = image->data;
    int width = image->width;
    int height = image->height;
    
    float expScale = -1.f / (d * d * 0.5f);
    
    float sin = sinf(angle * PI / 180.0f);
    float cos = cosf(angle * PI / 180.0f);
    
    /**360 to n*/
    float binPerRad = n / 360.f;
    
    float histWidth = SIFT_DESCRIPTOR_SCAE_FCTER * scale;
    
    int radius = (int) (roundf(histWidth * 1.4142135623730951f * (d + 1) * 0.5f));
    
    radius = min_value(radius, (int) sqrt(width * width + height * height));
    
    sin /= histWidth;
    cos /= histWidth;
    
    int len = (radius + radius + 1) * (radius + radius + 1);
    int histLen = (d + 2) * (d + 2) * (n + 2);
    
    float *buffer = (float *) malloc(sizeof(float) * (len * 7 + histLen));
    float *dx = buffer, *dy = dx + len, *mag = dy + len, *ori = mag + len, *weight = ori + len;
    float *rBin = weight + len, *cBin = rBin + len, *hist = cBin + len;
    
    memset(hist, 0, sizeof(float) * histLen);
    
    int realLen = 0;
    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
            float cRotate = cos * j - sin * i;
            float rRotate = sin * j + cos * i;
            
            float rBin0 = rRotate + d / 2.0f - 0.5f;
            float cBin0 = cRotate + d / 2.0f - 0.5f;
            
            int realX = x + j;
            int realY = y + i;
            
            if (rBin0 > -1 && rBin0 < d && cBin0 > -1 && cBin0 < d
                && realX > 0 && realX < width - 1 && realY > 0 && realY < height - 1) {
                dx[realLen] = (float) (data[realY * width + realX + 1] - data[realY * width + realX - 1]);
                dy[realLen] = (float) (data[(realY - 1) * width + realX] - data[(realY + 1) * width + realX]);
                
                rBin[realLen] = rBin0;
                cBin[realLen] = cBin0;
                
                weight[realLen] = (rRotate * rRotate + cRotate * cRotate) * expScale;
                
                realLen++;
            }
        }
    }
    
    for (int i = 0; i < realLen; ++i) {
        ori[i] = atan2f360(dy[i], dx[i]);
        mag[i] = sqrtf(dy[i] * dy[i] + dx[i] * dx[i]);
        weight[i] = exp(weight[i]);
    }
    
    for (int i = 0; i < realLen; ++i) {
        float rbin = rBin[i];
        float cbin = cBin[i];
        float obin = (ori[i] - angle) * binPerRad;
        
        float magnitude = mag[i] * weight[i];
        
        int r0 = (int) floorf(rbin);
        int c0 = (int) floorf(cbin);
        int o0 = (int) floorf(obin);
        
        rbin -= r0;
        cbin -= c0;
        obin -= o0;
        
        if (o0 < 0) {
            o0 += n;
        }
        if (o0 >= n) {
            o0 -= n;
        }
        
        float v_r1 = magnitude * rbin, v_r0 = magnitude - v_r1;
        float v_rc11 = v_r1 * cbin, v_rc10 = v_r1 - v_rc11;
        float v_rc01 = v_r0 * cbin, v_rc00 = v_r0 - v_rc01;
        float v_rco111 = v_rc11 * obin, v_rco110 = v_rc11 - v_rco111;
        float v_rco101 = v_rc10 * obin, v_rco100 = v_rc10 - v_rco101;
        float v_rco011 = v_rc01 * obin, v_rco010 = v_rc01 - v_rco011;
        float v_rco001 = v_rc00 * obin, v_rco000 = v_rc00 - v_rco001;
        
        int idx = ((r0 + 1) * (d + 2) + c0 + 1) * (n + 2) + o0;
        
        hist[idx] += v_rco000;
        hist[idx + 1] += v_rco001;
        hist[idx + (n + 2)] += v_rco010;
        hist[idx + (n + 3)] += v_rco011;
        hist[idx + (d + 2) * (n + 2)] += v_rco100;
        hist[idx + (d + 2) * (n + 2) + 1] += v_rco101;
        hist[idx + (d + 3) * (n + 2)] += v_rco110;
        hist[idx + (d + 3) * (n + 2) + 1] += v_rco111;
    }
    
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            int idx = ((i + 1) * (d + 2) + (j + 1)) * (n + 2);
            hist[idx] += hist[idx + n];
            hist[idx + 1] += hist[idx + n + 1];
            
            for (int k = 0; k < n; k++) {
                descriptor[(i * d + j) * n + k] = hist[idx + k];
            }
        }
    }
    
    len = d * d * n;
    float norm = 0;
    for (int i = 0; i < len; ++i) {
        norm += descriptor[i] * descriptor[i];
    }
    
    float threshold = sqrt(norm) * SIFT_DESCRIPTOR_MAGNITUDE_THRESHOLD;
    
    norm = 0;
    for (int i = 0; i < len; ++i) {
        descriptor[i] = min_value(threshold, descriptor[i]);
        norm += descriptor[i] * descriptor[i];
    }
    
    norm = SIFT_DESCRIPTOR_FCTOR / max_value(sqrt(norm), FLT_EPSILON);
    for (int i = 0; i < len; ++i) {
        descriptor[i] = (uint8_t) (descriptor[i] * norm);
    }
    
    free(buffer);
}

void calculateDescriptor(vector<vector<Mat<int16_t> *> > &gaussPyramid, vector<KeyPoint> &keyPoints, int d, int n) {
    for (auto iter = keyPoints.begin(); iter < keyPoints.end(); ++iter) {
        int octave = (*iter).octave;
        int octaveLayer = (*iter).octaveLayer;
        
        int x = (*iter).octaveX;
        int y = (*iter).octaveY;
        
        float angle = 360.f - (*iter).angle;
        float size = (*iter).octaveSize;
        
        (*iter).descriptor = (float *) malloc(sizeof(float) * d * d * n);
        
        calculateDescriptorOne(gaussPyramid[octave][octaveLayer],
                               x, y,
                               angle, size,
                               d, n,
                               (*iter).descriptor);
    }
}

/**
 * calculate the descriptor
 * @param gaussPyramid
 * @param keyPoints
 * @param d
 * @param n
 */
void calculateDescriptor1(vector<vector<Mat<int16_t> *> > &gaussPyramid, vector<KeyPoint> &keyPoints, int d, int n) {
    int size = keyPoints.size();
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex  = 0;
    
    /**remove the border*/
    int stride = max((int) (roundf(1.0f * size / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? size : min_value(start + stride, size);
        
        for (int i = start; i < end; ++i) {
            int octave = keyPoints[i].octave;
            int octaveLayer = keyPoints[i].octaveLayer;
            
            int x = keyPoints[i].octaveX;
            int y = keyPoints[i].octaveY;
            
            float angle = 360.f - keyPoints[i].angle;
            float size = keyPoints[i].octaveSize;
            
            keyPoints[i].descriptor = (float *) malloc(sizeof(float) * d * d * n);
            
            calculateDescriptorOne(gaussPyramid[octave][octaveLayer],
                                   x, y,
                                   angle, size,
                                   d, n,
                                   keyPoints[i].descriptor);
        }
    }
}

template<class T>
void deleteMatPyramid(vector<vector<SiftySifty::Mat<T> * >> &pyramid) {
    auto i1 = pyramid.begin();
    
    for (; i1 < pyramid.end(); ++i1) {
        auto i2 = (*i1).begin();
        
        for (; i2 < (*i1).end(); ++i2) {
            deleteMat<T>(*i2);
        }
        
        (*i1).clear();
    }
    
    pyramid.clear();
}

void sift(uint8_t *image, int width, int height,
          vector<SiftySifty::KeyPoint> &keyPoints,
          int octaveLayers,
          float sigma,
          float contrastThreshold,
          int edgeThreshold,
          bool doubleInitImage,
          int descriptorWidth,
          int descriptorHistBin) {
    
    Mat<int16_t> *base = initBaseImage(image, width, height, doubleInitImage, sigma, SIFT_IMAGE_SCALE_SHIFT);
   
    int octave = (int) (round(log(min(base->width, base->height)) / log(2) - 2));
    
    vector<vector<Mat<int16_t> *> > gaussPyramid = buildGaussPyramid(base, octave, octaveLayers, sigma);
    
    vector<vector<Mat<int16_t> *> > dogPyramid = buildDoGPyramid(gaussPyramid, octave, octaveLayers);
    
    findExtremaPoint(gaussPyramid,
                     dogPyramid,
                     keyPoints,
                     octave,
                     octaveLayers,
                     contrastThreshold,
                     edgeThreshold,
                     sigma);
    
    if (doubleInitImage) {
        resizeKeyPoints(keyPoints);
    }
    
    removeDoubleKeyPoints(keyPoints);

    /**calcaulate descip*/
    calculateDescriptor1(gaussPyramid, keyPoints, descriptorWidth, descriptorHistBin);
    
    deleteMatPyramid(gaussPyramid);
    deleteMatPyramid(dogPyramid);
}

void sift(uint8_t *image, int width, int height, vector<SiftySifty::KeyPoint> &keyPoints) {
    sift(image, width, height,
         keyPoints,
         SIFT_OCTAVE_LAYERS,
         SIFT_SIGMA,
         SIFT_CONTRAST_THRESHOLD,
         SIFT_EDGE_THESHOLD,
         SIFT_DOUBLE_INITED_IMAGE,
         SIFT_DESCRIPTOR_WIDTH,
         SIFT_DESCRIPTOR_HIST_BIN);
}

}












