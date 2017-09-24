/**
 * Created by yanyuanchi on 2017/1/16.
 */

#ifndef SIFT_UTILS_H
#define SIFT_UTILS_H

#include <stdlib.h>
#include <math.h>
#include <zconf.h>
#include <sys/time.h>
#include <string.h>

#include "structs.h"

#ifdef _OPENMP

#include <omp.h>

#endif

namespace SiftySifty {

#ifndef max_value
#define max_value(a, b) ((a) > ((b)) ? (a) : (b))
#endif

#ifndef min_value
#define min_value(a, b) (((a) < (b)) ? (a) : (b))
#endif

/**PI*/
static float PI = 3.1415926f;

static int HARDWARE_CPU_NUM = -1;

static int getHardwareCPUNum() {
    if (0 >= HARDWARE_CPU_NUM) {
        HARDWARE_CPU_NUM = static_cast<int>(sysconf(_SC_NPROCESSORS_CONF));
        
        if (0 >= HARDWARE_CPU_NUM) {
            HARDWARE_CPU_NUM = 4;
        }
    }
    
    return HARDWARE_CPU_NUM;
}

static long getCurrentTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

/**
 * get the angle of [0, 360)
 * @param y
 * @param x
 * @return
 */
static float atan2f360(float y, float x) {
    float aX = fabsf(x);
    float aY = fabsf(y);
    
    if (0 == x) {
        return (y > 0) ? 90 : 270;
    }
    
    if (0 == y) {
        return (x >= 0) ? 0 : 180;
    }
    
    float angle = atan2f(aY, aX) * 180.f / PI;
    
    if (x > 0) {
        return (y > 0) ? angle : (360 - angle);
    } else {
        return (y > 0) ? (180 - angle) : (180 + angle);
    }
}

template<class T>
Mat<T> *newMat(int width, int height) {
    if (0 >= width || 0 >= height) {
        return nullptr;
    }
    
    Mat<T> *mat = (Mat<T> *) malloc(sizeof(Mat<T>));
    mat->width = width;
    mat->height = height;
    mat->data = (T *) malloc(sizeof(T) * width * height);
    
    return mat;
}

template<class T>
void deleteMat(Mat<T> *mat) {
    if (nullptr != mat) {
        if (nullptr != mat->data) {
            free(mat->data);
        }
        
        free(mat);
    }
}

}
#endif //SIFT_UTILS_H
