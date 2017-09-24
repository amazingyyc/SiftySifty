/**
 * Created by yanyuanchi on 2017/4/8.
 */

#ifndef SIFTYSIFTY_IMAGEUTILS_H
#define SIFTYSIFTY_IMAGEUTILS_H

#include "structs.h"
#include "utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace SiftySifty {

/**
 * half sample the src Mat
 * @tparam T type
 * @param src src mat
 * @param dst dst mat
 * @return true/false
 */
template<class T>
void halfSampleMat(Mat<T> *src, Mat<T> *dst) {
    if (nullptr == src || nullptr == dst || nullptr == src->data || nullptr == dst->data) {
        return;
    }
    
    int srcWidth  = src->width;
    int srcHeight = src->height;
    
    int dstWidth  = dst->width;
    int dstHeight = dst->height;
    
    if ((srcWidth >> 1) != dstWidth || (srcHeight >> 1) != dstHeight) {
        return;
    }
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex  = 0;
    
    int stride = max_value((int) (roundf(1.0f * dstHeight / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end   = (threadIndex == (maxThreadNum - 1)) ? dstHeight : min_value(start + stride, dstHeight);
        
        T *srcData = src->data + (srcWidth * (start << 1));
        T *dstData = dst->data + (dstWidth * start);
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < dstWidth; ++x) {
                dstData[x] = srcData[(x << 1)];
            }
            
            srcData += (srcWidth + srcWidth);
            dstData += dstWidth;
        }
    }
}


/**
 * resize src to dst
 * @tparam T
 * @param src
 * @param dst
 * @return
 */
template<class T>
void resizeMat(Mat<T> *src, Mat<T> *dst) {
    if (nullptr == src || nullptr == dst || nullptr == src->data || nullptr == dst->data) {
        return;
    }
    
    int srcWidth = src->width;
    int srcHeight = src->height;
    
    int dstWidth = dst->width;
    int dstHeight = dst->height;
    
    T *srcData = src->data;
    T *dstData = dst->data;
    
    if (srcWidth == dstWidth && srcHeight == dstHeight) {
        memcpy(dstData, srcData, sizeof(T) * srcWidth * srcHeight);
        
        return;
    }
    
    int32_t shift = 22;
    int64_t scale = (1 << (shift >> 1));
    int64_t delta = (1 << (shift - 1));
    
    float xRatio = 1.f * (srcWidth - 1.f) / dstWidth;
    float yRatio = 1.f * (srcHeight - 1.f) / dstHeight;
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex  = 0;
    
    int stride = max_value((int) (roundf(1.0f * dstHeight / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? dstHeight : min_value(start + stride, dstHeight);
        
        T *dstOffsetData = dstData + start * dstWidth;
        
        for (int y = start; y < end; ++y) {
            float yOffset = (y + 0.5f) * yRatio;
            int yUp = (int) floorf(yOffset);
            
            yOffset -= yUp;
            
            int64_t multUp = (int64_t) (yOffset * scale);
            int64_t multDown = scale - multUp;
            
            for (int x = 0; x < dstWidth; ++x) {
                float xOffset = (x + 0.5f) * xRatio;
                int xLeft = (int) floorf(xOffset);
                
                xOffset -= xLeft;
                
                int64_t multLeft = (int64_t) (xOffset * scale);
                int64_t multRight = scale - multLeft;
                
                T *srcOffsetData = srcData + yUp * srcWidth + xLeft;
                
                dstOffsetData[x] = (T) ((srcOffsetData[0] * multRight * multDown
                                         + srcOffsetData[1] * multLeft * multDown
                                         + srcOffsetData[srcWidth] * multRight * multUp
                                         + srcOffsetData[srcWidth + 1] * multLeft * multUp
                                         + delta) >> shift);
            }
            
            dstOffsetData += dstWidth;
        }
    }
}

/**
 * resize the src to dst
 * @tparam T
 * @param src
 * @param srcWidth
 * @param srcHeight
 * @param dst
 * @param dstWidth
 * @param dstHeight
 * @return
 */
template<class T>
void resizeMat2(Mat<T> *srcMat, Mat<T> *dstMat) {
    if (nullptr == srcMat || nullptr == dstMat) {
        return;
    }
    
    T *src = srcMat->data;
    T *dst = dstMat->data;
    
    int srcWidth  = srcMat->width;
    int srcHeight = srcMat->height;
    int dstWidth  = dstMat->width;
    int dstHeight = dstMat->height;

    if (srcWidth == dstWidth && srcHeight == dstHeight) {
        memcpy(dst, src, sizeof(T) * srcWidth * srcHeight);
        return;
    }
    
    int32_t shift = 22;
    int64_t scale = (1 << (shift >> 1));
    int64_t delta = (1 << (shift - 1));
    
    /**
     * src = (dst + 0.5) * srcWidth / dstWidth - 0.5
     */
    float xRatio = 1.f * srcWidth / dstWidth;
    float yRatio = 1.f * srcHeight / dstHeight;
    
    int *xTable = (int*) malloc(sizeof(int) * 2 * dstWidth);
    int64_t *xMult = (int64_t*) malloc(sizeof(int64_t) * 2 * dstWidth);
    
    for (int x = 0; x < dstWidth; ++x) {
        float xOffset = (x + 0.5f) * xRatio - 0.5f;
        int xLeft;
        
        int64_t multLeft, multRight;
        
        if (0 >= xOffset) {
            xLeft = 0;
            multLeft = 0;
            multRight = scale;
        } else if (xOffset >= (srcWidth - 1)) {
            xLeft = srcWidth - 2;
            multLeft = scale;
            multRight = 0;
        } else {
            xLeft = (int) floorf(xOffset);
            
            xOffset -= xLeft;
            
            multLeft = (int64_t) (xOffset * scale);
            multRight = scale - multLeft;
        }
        
        xTable[(x << 1)] = xLeft;
        xTable[(x << 1) + 1] = xLeft + 1;
        xMult[(x << 1)] = multLeft;
        xMult[(x << 1) + 1] = multRight;
    }
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex  = 0;
    
    int stride = max_value((int) (roundf(1.0f * dstHeight / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end   = (threadIndex == (maxThreadNum - 1)) ? dstHeight : min_value(start + stride, dstHeight);
        
        T *dstData = dst + start * dstWidth;
        
        for (int y = start; y < end; ++y) {
            float yOffset = (y + 0.5f) * yRatio - 0.5f;
            int yUp;
            int64_t multUp, multDown;
            
            if (0 >= yOffset) {
                yUp = 0;
                multUp = 0;
                multDown = scale;
            } else if (yOffset >= (srcHeight - 1)) {
                yUp = srcHeight - 2;
                multUp = scale;
                multDown = 0;
            } else {
                yUp = (int) floorf(yOffset);
                yOffset -= yUp;
                
                multUp = (int64_t) (yOffset * scale);
                multDown = scale - multUp;
            }
            
            T *upSrc = src + yUp * srcWidth;
            T *downSrc = upSrc + srcWidth;
            
            for (int x = 0; x < dstWidth; ++x) {
                int x2 = (x << 1);
                
                dstData[x] = (T)((((upSrc[xTable[x2]] * xMult[x2+1]
                                    + upSrc[xTable[x2+1]]*xMult[x2]) * multDown
                                   + (downSrc[xTable[x2]] * xMult[x2+1]
                                      + downSrc[xTable[x2+1]]*xMult[x2]) * multUp)
                                  + delta) >> shift);
            }
            
            dstData += dstWidth;
        }
    }
    
    free(xTable);
    free(xMult);
}

template<class T1, class T2>
void scaleMatByScale(T1 *src, T2 *dst, int width, int height, int scale) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex  = 0;
    
    int stride = max_value((int) (roundf(1.0f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : min_value(start + stride, height);
        
        T1 *srcData = src + start * width;
        T2 *dstData = dst + start * width;
        
        int64_t length = (end - start) * width;
        int64_t limit = length - 3;
        
        int64_t i = 0;
        for (; i < limit; i += 4) {
            dstData[i] = (srcData[i] * scale);
            dstData[i + 1] = (srcData[i + 1] * scale);
            dstData[i + 2] = (srcData[i + 2] * scale);
            dstData[i + 3] = (srcData[i + 3] * scale);
        }
        
        for (; i < length; ++i) {
            dstData[i] = (srcData[i] * scale);
        }
    }
}

template<class T1, class T2>
void scaleMatByShift(T1 *src, T2 *dst, int width, int height, int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = max_value((int) (roundf(1.0f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : min_value(start + stride, height);
        
        T1 *srcData = src + start * width;
        T2 *dstData = dst + start * width;
        
        int64_t length = (end - start) * width;
        int64_t limit = length - 3;
        
        int64_t i = 0;
        for (; i < limit; i += 4) {
            dstData[i] = (srcData[i] << shift);
            dstData[i + 1] = (srcData[i + 1] << shift);
            dstData[i + 2] = (srcData[i + 2] << shift);
            dstData[i + 3] = (srcData[i + 3] << shift);
        }
        
        for (; i < length; ++i) {
            dstData[i] = (srcData[i] << shift);
        }
    }
}

template<class T>
void subMat(T *src1, T *src2, T *dst, int width, int height) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex  = 0;
    
    int stride = max_value((int) (roundf(1.0f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : min_value(start + stride, height);
        
        T *src1Data = src1 + start * width;
        T *src2Data = src2 + start * width;
        T *dstData  = dst + start * width;
        
        int64_t length = (end - start) * width;
        int64_t limit  = length - 3;
        
        int64_t i = 0;
        for (; i < limit; i += 4) {
            dstData[i] = src1Data[i] - src2Data[i];
            dstData[i + 1] = src1Data[i + 1] - src2Data[i + 1];
            dstData[i + 2] = src1Data[i + 2] - src2Data[i + 2];
            dstData[i + 3] = src1Data[i + 3] - src2Data[i + 3];
        }
        
        for (; i < length; ++i) {
            dstData[i] = src1Data[i] - src2Data[i];
        }
    }
}

}


#endif //SIFTYSIFTY_IMAGEUTILS_H
