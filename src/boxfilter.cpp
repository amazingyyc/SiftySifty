/**
 * Created by yanyuanchi on 2017/3/21.
 */
#include <math.h>

#include "utils.h"
#include "filter.h"
#include "boxfilter.h"

namespace SiftySifty {

/**
 * boxFilter in row
 * dst[i] = sum(src[i - radius, i + radius]) / (2 * radius + 1)
 * scale = (int) (1.0f / (2 * radius + 1) * (1 << shift))
 * delta = (1 << (shift - 1))
 *
 * dst[i] = (sum(src[i - radius, i + radius]) * scale + delta) >> shift
 *
 * the src have border with radius's cols in left and radius's cols in right
 *
 * T can be uint8_t/int8_t/uint16_t/int16_t
 */
template<class T>
void boxFilterRow(T *src, T *dst, int width, int height, int radius, int scale, int delta, int shift) {
    int radius2 = radius + radius;
    int size = radius2 + 1;
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = max_value((int) (roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : min_value(start + stride, height);
        
        T *srcData = src + start * (width + radius2);
        T *dstData = dst + start * width;
        
        int64_t sum;
        
        for (int y = start; y < end; ++y) {
            sum = 0;
            
            for (int x = 0; x < size; ++x) {
                sum += srcData[x];
            }
            
            dstData[0] = (T) ((sum * scale + delta) >> shift);
            
            for (int x = 1; x < width; ++x) {
                sum += srcData[x + radius2] - srcData[x - 1];
                
                dstData[x] = (T) ((sum * scale + delta) >> shift);
            }
            
            srcData += (width + radius2);
            dstData += width;
        }
    }
}

/**
 * the src with radius rows int top and bottom
 */
template<class T>
void boxFilterCol(T *src, T *dst, int width, int height, int radius, int scale, int delta, int shift) {
    int radius2 = radius + radius;
    int size = radius2 + 1;
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    /**split the width's cols to maxThreadNum's thread*/
    int stride = max_value((int) (roundf(1.f * width / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? width : min_value(width, start + stride);
        
        int range = end - start;
        
        int interval = radius2 * width;
        
        int64_t *sum = (int64_t *) malloc(sizeof(int64_t) * range);
        memset(sum, 0, sizeof(int64_t) * range);
        
        T *srcData = src + start;
        T *dstData = dst + start;
        
        for (int y = 0; y < radius2; ++y) {
            for (int x = 0; x < range; ++x) {
                sum[x] += srcData[x];
            }
            
            srcData += width;
        }
        
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < range; ++x) {
                sum[x] += srcData[x];
                
                dstData[x] = (T) ((sum[x] * scale + delta) >> shift);
                
                sum[x] -= srcData[x - interval];
            }
            
            srcData += width;
            dstData += width;
        }
        
        free(sum);
    }
}

/**
 * box blur
 * first blur int row
 * than blur int cols
 */
template<class T>
void boxFilter(T *src, T *dst, int width, int height, int radius) {
    if (0 >= radius) {
        memcpy(dst, src, sizeof(T) * width * height);
        
        return;
    }
    
    int radius2 = radius + radius;
    int size = radius2 + 1;
    
    int scale = (int) (1.0 / size * FILTER_SCALE);
    
    T *srcTemp = (T *) malloc(sizeof(T) * (width + radius2) * height);
    
    T *srcData = src;
    T *srcTempData = srcTemp;
    
    /**copy memory to srcMediate*/
    for (int y = 0; y < height; ++y) {
        std::fill(srcTempData, srcTempData + radius, srcData[0]);
        
        memcpy(srcTempData + radius, srcData, sizeof(T) * width);
        
        std::fill(srcTempData + radius + width,
                  srcTempData + radius2 + width,
                  srcData[width - 1]);
        
        srcData += width;
        srcTempData += (width + radius2);
    }
    
    T *dstTemp = (T *) malloc(sizeof(T) * width * (height + radius2));
    
    /**blur in row*/
    boxFilterRow<T>(srcTemp, dstTemp + (radius * width),
                    width, height, radius,
                    scale, FILTER_DELTA, FILTER_SHIFT);
    
    for (int y = 0; y < radius; ++y) {
        memcpy(dstTemp + y * width, dstTemp + radius * width, sizeof(T) * width);
        memcpy(dstTemp + (radius + height + y) * width, dstTemp + (radius + height - 1) * width, sizeof(T) * width);
    }
    
    boxFilterCol<T>(dstTemp, dst, width, height, radius, scale, FILTER_DELTA, FILTER_SHIFT);
    
    free(srcTemp);
    free(dstTemp);
}

/**
 * box blur
 */
void boxFilter(int16_t *src, int16_t *dst, int width, int height, int radius) {
    boxFilter<int16_t>(src, dst, width, height, radius);
}

}




















