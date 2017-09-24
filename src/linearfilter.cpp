/**
 * Created by yanyuanchi on 2017/1/18.
 */

#include <cmath>

#include "utils.h"
#include "linearfilter.h"

namespace SiftySifty {
void linearFilterHorizonByKernel7(uint8_t *src,
                                  uint8_t *dst,
                                  int width,
                                  int height,
                                  int (*mult)[256],
                                  int delta,
                                  int shift) {
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * (width + 6);
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                sum += mult[1][realSrcData[1]];
                sum += mult[2][realSrcData[2]];
                sum += mult[3][realSrcData[3]];
                sum += mult[4][realSrcData[4]];
                sum += mult[5][realSrcData[5]];
                sum += mult[6][realSrcData[6]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += (width + 6);
            dstData += width;
        }
    }
}

void linearFilterHorizonByKernel9(uint8_t *src,
                                  uint8_t *dst,
                                  int width,
                                  int height,
                                  int (*mult)[256],
                                  int delta,
                                  int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * (width + 8);
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                sum += mult[1][realSrcData[1]];
                sum += mult[2][realSrcData[2]];
                sum += mult[3][realSrcData[3]];
                sum += mult[4][realSrcData[4]];
                sum += mult[5][realSrcData[5]];
                sum += mult[6][realSrcData[6]];
                sum += mult[7][realSrcData[7]];
                sum += mult[8][realSrcData[8]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += (width + 8);
            dstData += width;
        }
    }
}

void linearFilterHorizonByKernel11(uint8_t *src,
                                   uint8_t *dst,
                                   int width,
                                   int height,
                                   int (*mult)[256],
                                   int delta,
                                   int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * (width + 10);
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                sum += mult[1][realSrcData[1]];
                sum += mult[2][realSrcData[2]];
                sum += mult[3][realSrcData[3]];
                sum += mult[4][realSrcData[4]];
                sum += mult[5][realSrcData[5]];
                sum += mult[6][realSrcData[6]];
                sum += mult[7][realSrcData[7]];
                sum += mult[8][realSrcData[8]];
                sum += mult[9][realSrcData[9]];
                sum += mult[10][realSrcData[10]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += (width + 10);
            dstData += width;
        }
    }
}


void linearFilterHorizonByKernel13(uint8_t *src,
                                   uint8_t *dst,
                                   int width,
                                   int height,
                                   int (*mult)[256],
                                   int delta,
                                   int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * (width + 12);
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                sum += mult[1][realSrcData[1]];
                sum += mult[2][realSrcData[2]];
                sum += mult[3][realSrcData[3]];
                sum += mult[4][realSrcData[4]];
                sum += mult[5][realSrcData[5]];
                sum += mult[6][realSrcData[6]];
                sum += mult[7][realSrcData[7]];
                sum += mult[8][realSrcData[8]];
                sum += mult[9][realSrcData[9]];
                sum += mult[10][realSrcData[10]];
                sum += mult[11][realSrcData[11]];
                sum += mult[12][realSrcData[12]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += (width + 12);
            dstData += width;
        }
    }
}

void linearFilterHorizonByKernel15(uint8_t *src,
                                   uint8_t *dst,
                                   int width,
                                   int height,
                                   int (*mult)[256],
                                   int delta,
                                   int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * (width + 14);
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                sum += mult[1][realSrcData[1]];
                sum += mult[2][realSrcData[2]];
                sum += mult[3][realSrcData[3]];
                sum += mult[4][realSrcData[4]];
                sum += mult[5][realSrcData[5]];
                sum += mult[6][realSrcData[6]];
                sum += mult[7][realSrcData[7]];
                sum += mult[8][realSrcData[8]];
                sum += mult[9][realSrcData[9]];
                sum += mult[10][realSrcData[10]];
                sum += mult[11][realSrcData[11]];
                sum += mult[12][realSrcData[12]];
                sum += mult[13][realSrcData[13]];
                sum += mult[14][realSrcData[14]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += (width + 14);
            dstData += width;
        }
    }
}

void linearFilterHorizonByKernel17(uint8_t *src,
                                   uint8_t *dst,
                                   int width,
                                   int height,
                                   int (*mult)[256],
                                   int delta,
                                   int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * (width + 16);
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                sum += mult[1][realSrcData[1]];
                sum += mult[2][realSrcData[2]];
                sum += mult[3][realSrcData[3]];
                sum += mult[4][realSrcData[4]];
                sum += mult[5][realSrcData[5]];
                sum += mult[6][realSrcData[6]];
                sum += mult[7][realSrcData[7]];
                sum += mult[8][realSrcData[8]];
                sum += mult[9][realSrcData[9]];
                sum += mult[10][realSrcData[10]];
                sum += mult[11][realSrcData[11]];
                sum += mult[12][realSrcData[12]];
                sum += mult[13][realSrcData[13]];
                sum += mult[14][realSrcData[14]];
                sum += mult[15][realSrcData[15]];
                sum += mult[16][realSrcData[16]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += (width + 16);
            dstData += width;
        }
    }
}

void linearFilterHorizonByKernel19(uint8_t *src,
                                   uint8_t *dst,
                                   int width,
                                   int height,
                                   int (*mult)[256],
                                   int delta,
                                   int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * (width + 18);
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                sum += mult[1][realSrcData[1]];
                sum += mult[2][realSrcData[2]];
                sum += mult[3][realSrcData[3]];
                sum += mult[4][realSrcData[4]];
                sum += mult[5][realSrcData[5]];
                sum += mult[6][realSrcData[6]];
                sum += mult[7][realSrcData[7]];
                sum += mult[8][realSrcData[8]];
                sum += mult[9][realSrcData[9]];
                sum += mult[10][realSrcData[10]];
                sum += mult[11][realSrcData[11]];
                sum += mult[12][realSrcData[12]];
                sum += mult[13][realSrcData[13]];
                sum += mult[14][realSrcData[14]];
                sum += mult[15][realSrcData[15]];
                sum += mult[16][realSrcData[16]];
                sum += mult[17][realSrcData[17]];
                sum += mult[18][realSrcData[18]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += (width + 18);
            dstData += width;
        }
    }
}

void linearFilterHorizonByKernel(uint8_t *src,
                                 uint8_t *dst,
                                 int width,
                                 int height,
                                 int (*mult)[256],
                                 int delta,
                                 int shift,
                                 int size) {
    int radius = (size - 1) / 2;
    
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * (width + radius + radius);
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                for (int k = 0; k < size; ++k) {
                    sum += mult[k][realSrcData[k]];
                }
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += width + radius + radius;
            dstData += width;
        }
    }
}

void linearFilterHorizon(uint8_t *src,
                         uint8_t *dst,
                         int width,
                         int height,
                         int (*mult)[256],
                         int delta,
                         int shift,
                         int size) {
    int radius = (size - 1) / 2;
    
    uint8_t *tmp = new uint8_t[(width + radius + radius) * height];
    
    /**
     * --------------------------------
     * |       |              |       |
     * |radius |              |radius |
     * |       |              |       |
     * |       |              |       |
     * |       |              |       |
     * |       |              |       |
     * --------------------------------
     */
    uint8_t *srcOffset = src;
    uint8_t *tmpOffset = tmp;
    
    for (int y = 0; y < height; ++y) {
        std::fill(tmpOffset,
                  tmpOffset + radius,
                  srcOffset[0]);
        
        memcpy(tmpOffset + radius,
               srcOffset,
               sizeof(uint8_t) * width);
        
        std::fill(tmpOffset + radius + width,
                  tmpOffset + +radius + radius + width,
                  srcOffset[width - 1]);
        
        srcOffset += width;
        tmpOffset += width + radius + radius;
    }
    
    if (7 == size) {
        linearFilterHorizonByKernel7(tmp, dst, width, height, mult, delta, shift);
    } else if (9 == size) {
        linearFilterHorizonByKernel9(tmp, dst, width, height, mult, delta, shift);
    } else if (11 == size) {
        linearFilterHorizonByKernel11(tmp, dst, width, height, mult, delta, shift);
    } else if (13 == size) {
        linearFilterHorizonByKernel13(tmp, dst, width, height, mult, delta, shift);
    } else if (15 == size) {
        linearFilterHorizonByKernel15(tmp, dst, width, height, mult, delta, shift);
    } else if (17 == size) {
        linearFilterHorizonByKernel17(tmp, dst, width, height, mult, delta, shift);
    } else if (19 == size) {
        linearFilterHorizonByKernel19(tmp, dst, width, height, mult, delta, shift);
    } else {
        linearFilterHorizonByKernel(tmp, dst, width, height, mult, delta, shift, size);
    }
    
    delete[] tmp;
}


void linearFilterVerticalByKernel7(uint8_t *src,
                                   uint8_t *dst,
                                   int width,
                                   int height,
                                   int (*mult)[256],
                                   int delta,
                                   int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * width;
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[1][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[2][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[3][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[4][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[5][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[6][realSrcData[0]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += width;
            dstData += width;
        }
    }
}

void linearFilterVerticalByKernel9(uint8_t *src,
                                   uint8_t *dst,
                                   int width,
                                   int height,
                                   int (*mult)[256],
                                   int delta,
                                   int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * width;
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[1][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[2][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[3][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[4][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[5][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[6][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[7][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[8][realSrcData[0]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += width;
            dstData += width;
        }
    }
}

void linearFilterVerticalByKernel11(uint8_t *src,
                                    uint8_t *dst,
                                    int width,
                                    int height,
                                    int (*mult)[256],
                                    int delta,
                                    int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * width;
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[1][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[2][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[3][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[4][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[5][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[6][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[7][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[8][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[9][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[10][realSrcData[0]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += width;
            dstData += width;
        }
    }
}

void linearFilterVerticalByKernel13(uint8_t *src,
                                    uint8_t *dst,
                                    int width,
                                    int height,
                                    int (*mult)[256],
                                    int delta,
                                    int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * width;
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[1][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[2][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[3][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[4][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[5][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[6][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[7][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[8][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[9][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[10][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[11][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[12][realSrcData[0]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += width;
            dstData += width;
        }
    }
}

void linearFilterVerticalByKernel15(uint8_t *src,
                                    uint8_t *dst,
                                    int width,
                                    int height,
                                    int (*mult)[256],
                                    int delta,
                                    int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * width;
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[1][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[2][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[3][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[4][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[5][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[6][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[7][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[8][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[9][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[10][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[11][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[12][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[13][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[14][realSrcData[0]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += width;
            dstData += width;
        }
    }
}

void linearFilterVerticalByKernel17(uint8_t *src,
                                    uint8_t *dst,
                                    int width,
                                    int height,
                                    int (*mult)[256],
                                    int delta,
                                    int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * width;
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[1][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[2][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[3][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[4][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[5][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[6][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[7][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[8][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[9][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[10][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[11][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[12][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[13][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[14][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[15][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[16][realSrcData[0]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += width;
            dstData += width;
        }
    }
}

void linearFilterVerticalByKernel19(uint8_t *src,
                                    uint8_t *dst,
                                    int width,
                                    int height,
                                    int (*mult)[256],
                                    int delta,
                                    int shift) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * width;
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int sum;
        
        for (int y = start; y < end; ++y) {
            for (int x = 0; x < width; ++x) {
                realSrcData = srcData + x;
                sum = delta;
                
                sum += mult[0][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[1][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[2][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[3][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[4][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[5][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[6][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[7][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[8][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[9][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[10][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[11][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[12][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[13][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[14][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[15][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[16][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[17][realSrcData[0]];
                
                realSrcData += width;
                sum += mult[18][realSrcData[0]];
                
                dstData[x] = static_cast<uint8_t>((sum + delta) >> shift);
            }
            
            srcData += width;
            dstData += width;
        }
    }
}

void linearFilterVerticalByKernel(uint8_t *src,
                                  uint8_t *dst,
                                  int width,
                                  int height,
                                  int (*mult)[256],
                                  int delta,
                                  int shift,
                                  int size) {
    int maxThreadNum = getHardwareCPUNum();
    int threadIndex = 0;
    
    int stride = std::max(static_cast<int>(roundf(1.f * height / maxThreadNum)), 1);

#pragma omp parallel for private(threadIndex)
    for (threadIndex = 0; threadIndex < maxThreadNum; ++threadIndex) {
        int start = threadIndex * stride;
        int end = (threadIndex == (maxThreadNum - 1)) ? height : std::min(start + stride, height);
        
        uint8_t *srcData = src + start * width;
        uint8_t *dstData = dst + start * width;
        
        uint8_t *realSrcData;
        
        int *sum = (int *) malloc(sizeof(int) * width);
        
        for (int y = start; y < end; ++y) {
            std::fill(sum, sum + width, delta);
            
            realSrcData = srcData;
            
            for (int k = 0; k < size; ++k) {
                for (int x = 0; x < width; ++x) {
                    sum[x] += mult[k][realSrcData[x]];
                }
                
                realSrcData += width;
            }
            
            for (int x = 0; x < width; ++x) {
                dstData[x] = static_cast<uint8_t>(((sum[x] + delta) >> shift));
            }
            
            srcData += width;
            dstData += width;
        }
        
        free(sum);
    }
}

void linearFilterVertical(uint8_t *src,
                          uint8_t *dst,
                          int width,
                          int height,
                          int (*mult)[256],
                          int delta,
                          int shift,
                          int size) {
    int radius = (size - 1) / 2;
    uint8_t *tmp = new uint8_t[width * (height + radius + radius)];
    
    memcpy(tmp + radius * width, src, sizeof(uint8_t) * width * height);
    
    for (int i = 0; i < radius; ++i) {
        memcpy(tmp + i * width, src, sizeof(uint8_t) * width);
        memcpy(tmp + (radius + height + i) * width, src + (height - 1) * width, sizeof(uint8_t) * width);
    }
    
    if (7 == size) {
        linearFilterVerticalByKernel7(tmp, dst, width, height, mult, delta, shift);
    } else if (9 == size) {
        linearFilterVerticalByKernel9(tmp, dst, width, height, mult, delta, shift);
    } else if (11 == size) {
        linearFilterVerticalByKernel11(tmp, dst, width, height, mult, delta, shift);
    } else if (13 == size) {
        linearFilterVerticalByKernel13(tmp, dst, width, height, mult, delta, shift);
    } else if (15 == size) {
        linearFilterVerticalByKernel15(tmp, dst, width, height, mult, delta, shift);
    } else if (17 == size) {
        linearFilterVerticalByKernel17(tmp, dst, width, height, mult, delta, shift);
    } else if (19 == size) {
        linearFilterVerticalByKernel19(tmp, dst, width, height, mult, delta, shift);
    } else {
        linearFilterVerticalByKernel(tmp, dst, width, height, mult, delta, shift, size);
    }
    
    delete[] tmp;
}
    
}



























