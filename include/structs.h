/**
 * Created by yanyuanchi on 2017/1/16.
 */

#ifndef SIFT_STRUCTS_H_H
#define SIFT_STRUCTS_H_H

#include <iostream>
#include <cstdint>

namespace SiftySifty {

typedef struct KeyPoint {
    /**the coordinate of keyPoint*/
    float x;
    float y;
    
    /**the coordinate in pyramid*/
    int octaveX;
    int octaveY;
    
    /**nothing*/
    float score;
    
    float response;
    
    int octave;
    
    int octaveLayer;
    
    float octaveLayersShift;
    
    /**the scale of the keyPoint*/
    float size;
    
    /**the scale of the keyPoint, in the pyramid*/
    float octaveSize;
    
    /**the direction of the keypoint*/
    float angle;
    
    /**the descriptor of the keyPoint*/
    float *descriptor;
} KeyPoint;

template<class T>
struct Mat {
    T *data;
    
    int width;
    int height;
};

}
#endif //SIFT_STRUCTS_H_H
