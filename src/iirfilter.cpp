/**
 * Created by yanyuanchi on 2017/3/23.
 */

#include "utils.h"
#include "filter.h"
#include "iirfilter.h"

#ifdef _OPENMP

#include <omp.h>

#endif

namespace SiftySifty {
/**
 * use iir to filter the src
 * ref:"Recursive Implementation of the gaussian filter."
 * w[n] = (B * input[n] + b1 * w[n-1] + b2 * w[n-2] + b3 * w[n-3] + delta) >> shift
 */
template<class T>
void IIRFilterRow(T *src, T *dst,
                  const int width, const int height,
                  const int32_t B,
                  const int32_t b1, const int32_t b2, const int32_t b3,
                  const int32_t delta, const int32_t shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = max_value((int) (roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : min_value(start + stride, height);
        
        T *srcData = src + start * width;
        T *dstData = dst + start * width;
        
        int size;
        
        T *w = (T *) malloc(sizeof(T) * (width + 3));
        T tail;
        
        for (int y = start; y < end; ++y) {
            size = width - 1;
            
            w[0] = w[1] = w[2] = srcData[0];
            
            for (int x = 0, n = 3; x <= size; ++x, ++n) {
                w[n] = (T) ((B * srcData[x] + b1 * w[n - 1] + b2 * w[n - 2] + b3 * w[n - 3] + delta) >> shift);
            }
            
            tail = w[size + 3];
            
            dstData[size] = (T) ((B * w[size + 3] + b1 * tail + b2 * tail + b3 * tail + delta) >> shift);
            size--;
            dstData[size] = (T) ((B * w[size + 3] + b1 * dstData[size + 1] + b2 * tail + b3 * tail + delta) >> shift);
            size--;
            dstData[size] = (T) ((B * w[size + 3] + b1 * dstData[size + 1] + b2 * dstData[size + 2] + b3 * tail + delta)
                    >> shift);
            size--;
            
            for (int x = size; x >= 0; --x) {
                dstData[x] = (T) (
                        (B * w[x + 3] + b1 * dstData[x + 1] + b2 * dstData[x + 2] + b3 * dstData[x + 3] + delta)
                                >> shift);
            }
            
            srcData += width;
            dstData += width;
        }
        
        free(w);
    }
}

/**
 * filter on col same as the row
 * @tparam T
 * @param src
 * @param dst
 * @param width
 * @param height
 * @param delta
 * @param shift
 * @param B
 * @param b1
 * @param b2
 * @param b3
 */
template<class T>
void IIRFilterCol(T *src, T *dst,
                  const int width, const int height,
                  const int32_t B,
                  const int32_t b1, const int32_t b2, const int32_t b3,
                  const int32_t delta, const int32_t shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = max_value((int) (roundf(1.f * width / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? width : min_value(start + stride, width);
        int range = end - start;
        
        T *w = (T *) malloc(sizeof(T) * range * (height + 3));
        
        T *srcData = src + start;
        T *dstData = dst + start;
        T *wData = w + 3 * range;
        
        T *srcOffsetData;
        T *dstOffsetData;
        T *wOffsetData;
        
        memcpy(w, srcData, sizeof(T) * range);
        memcpy(w + range, srcData, sizeof(T) * range);
        memcpy(w + 2 * range, srcData, sizeof(T) * range);
        
        int off1 = -range;
        int off2 = off1 - range;
        int off3 = off2 - range;
        
        int size = height - 1;
        
        /**forward pass*/
        for (int y = 0; y <= size; ++y) {
            srcOffsetData = srcData;
            wOffsetData = wData;
            
            for (int x = 0; x < range; ++x) {
                wOffsetData[0] = (T) ((B * srcOffsetData[0] + b1 * wOffsetData[off1] + b2 * wOffsetData[off2] +
                                       b3 * wOffsetData[off3] + delta) >> shift);
                
                srcOffsetData++;
                wOffsetData++;
            }
            
            srcData += width;
            wData += range;
        }
        
        /**backward pass*/
        T *tail = (T *) malloc(sizeof(T) * range);
        memcpy(tail, w + range * (size + 3), sizeof(T) * range);
        
        off1 = width;
        off2 = off1 + width;
        off3 = off2 + width;
        
        dstData = dst + start + size * width;
        wData = w + (size + 3) * range;
        
        dstOffsetData = dstData;
        wOffsetData = wData;
        
        for (int x = 0; x < range; ++x) {
            dstOffsetData[0] = (T) ((B * wOffsetData[0] + b1 * tail[x] + b2 * tail[x] + b3 * tail[x] + delta) >> shift);
            
            dstOffsetData++;
            wOffsetData++;
        }
        
        dstData -= width;
        wData -= range;
        
        dstOffsetData = dstData;
        wOffsetData = wData;
        
        for (int x = 0; x < range; ++x) {
            dstOffsetData[0] = (T) ((B * wOffsetData[0] + b1 * dstOffsetData[off1] + b2 * tail[x] + b3 * tail[x]
                                     + delta) >> shift);
            
            dstOffsetData++;
            wOffsetData++;
        }
        
        dstData -= width;
        wData -= range;
        
        dstOffsetData = dstData;
        wOffsetData = wData;
        
        for (int x = 0; x < range; ++x) {
            dstOffsetData[0] = (T) (
                    (B * wOffsetData[0] + b1 * dstOffsetData[off1] + b2 * dstOffsetData[off2] + b3 * tail[x]
                     + delta) >> shift);
            
            dstOffsetData++;
            wOffsetData++;
        }
        
        dstData -= width;
        wData -= range;
        
        for (int y = size - 3; y >= 0; --y) {
            dstOffsetData = dstData;
            wOffsetData = wData;
            
            for (int x = 0; x < range; ++x) {
                dstOffsetData[0] = (T) ((B * wOffsetData[0] + b1 * dstOffsetData[off1] + b2 * dstOffsetData[off2]
                                         + b3 * dstOffsetData[off3] + delta) >> shift);
                
                dstOffsetData++;
                wOffsetData++;
            }
            
            dstData -= width;
            wData -= range;
        }
        
        free(tail);
        free(w);
    }
}

template<class T>
void IIRFilter(T *src, T *dst, const int width, const int height, const float sigma) {
    if (nullptr == src || nullptr == dst || 3 > width || 3 > height || 0 > sigma) {
        return;
    }
    
    double_t q, q2, q3;
    
    if (sigma >= 2.5) {
        q = 0.98711 * sigma - 0.96330;
    } else if (sigma >= 0.5 && sigma < 2.5) {
        q = 3.97156 - 4.14554 * sqrt(1.0 - 0.26891 * sigma);
    } else {
        q = 0.1147705018520355224609375;
    }
    
    q2 = q * q;
    q3 = q * q2;
    
    double_t db0 = 1.57825 + 2.44413 * q + 1.4281 * q2 + 0.422205 * q3;
    double_t db1 = 2.44413 * q + 2.85619 * q2 + 1.26661 * q3;
    double_t db2 = -(1.4281 * q2 + 1.26661 * q3);
    double_t db3 = 0.4222205 * q3;
    
    double_t dB = 1.0 - (db1 + db2 + db3) / db0;
    
    int32_t B  = (int32_t) (dB * FILTER_SCALE);
    int32_t b1 = (int32_t) (db1 / db0 * FILTER_SCALE);
    int32_t b2 = (int32_t) (db2 / db0 * FILTER_SCALE);
    int32_t b3 = (int32_t) (db3 / db0 * FILTER_SCALE);
    
    T *tmp = (T *) malloc(sizeof(T) * width * height);
    
    IIRFilterRow(src, tmp, width, height, B, b1, b2, b3, FILTER_DELTA, FILTER_SHIFT);
    IIRFilterCol(tmp, dst, width, height, B, b1, b2, b3, FILTER_DELTA, FILTER_SHIFT);
    
    free(tmp);
}

void IIRFilter(int16_t *src, int16_t *dst, int width, int height, float sigma) {
    IIRFilter<int16_t>(src, dst, width, height, sigma);
}

void IIRFilter(uint8_t *src, uint8_t *dst, int width, int height, float sigma) {
    IIRFilter<uint8_t>(src, dst, width, height, sigma);
}

}







