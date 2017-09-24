#include <iostream>

#include "siftysiftytest.h"

using namespace std;
using namespace SiftySifty;

#define IMAGE_PATH "../data/zgr.jpg"
#define IMAGE_PATH2 "../data/zgr2.jpg"
#define IMAGE_PATH_BIG "../data/zgrbig.jpg"

int main() {
    if (false) {
        drawKeyPoint(IMAGE_PATH);
    }
    
    if (false) {
        drawKeyPointCmpToOpenCV(IMAGE_PATH);
    }
    
    if (false) {
        matchKeyPoint(IMAGE_PATH, IMAGE_PATH);
    }
    
    if (true) {
        matchKeyPoint(IMAGE_PATH, IMAGE_PATH2);
    }
    
    if (false) {
        matchKeyPointSiftySiftyWithOpenCV(IMAGE_PATH, IMAGE_PATH);
    }
    
    if (false) {
        testSpeedSiftySiftyAndOpenCV(IMAGE_PATH_BIG);
    }
    
    return 0;
}









