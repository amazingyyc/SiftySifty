/**
 * Created by yanyuanchi on 2017/4/3.
 */
#ifndef SIFTYSIFTY_SIFTYTEST_H
#define SIFTYSIFTY_SIFTYTEST_H

#include <iostream>

using namespace std;

namespace SiftySifty {

/**
 * use the opencv to draw the siftysifty keypoint
 * @param path
 */
void drawKeyPoint(string path);

/**
 * draw the keypoint on the same pic by SiftySifty and OpenCV
 * @param path
 */
void drawKeyPointCmpToOpenCV(string path);

/**
 * math the keypoint that extracted by SiftySifty
 * @param path
 */
void matchKeyPoint(string path1, string path2);

/**
 * match the siftysifty keypoint with opencv keypoint
 * @param path
 */
void matchKeyPointSiftySiftyWithOpenCV(string path1, string path2);

/**
 * test the speed of SiftySifty and OpenCV
 * @param path
 */
void testSpeedSiftySiftyAndOpenCV(string path);

}

#endif //SIFTYSIFTY_SIFTYTEST_H








