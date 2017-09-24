/**
 * Created by yanyuanchi on 2017/4/3.
 */
#include <iostream>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <opencv2/imgproc.hpp>

#include "utils.h"
#include "siftysifty.h"
#include "siftysiftytest.h"

using namespace std;
using namespace cv;

namespace SiftySifty {
/**
 * use the opencv to draw the siftysifty keypoint
 * @param path
 */
void drawKeyPoint(string path)
{
    /**read the pic*/
    cv::Mat originImage = imread(path);
    cv::Mat grayImage;
    cvtColor(originImage, grayImage, CV_RGB2GRAY);

    vector<SiftySifty::KeyPoint> keyPoints;
    SiftySifty::sift(grayImage.data, grayImage.cols, grayImage.rows, keyPoints);

    vector<cv::KeyPoint> opencvKeyPoints(keyPoints.size());
    for (int i = 0; i < keyPoints.size(); ++i) {
        opencvKeyPoints[i].pt.x = keyPoints[i].x;
        opencvKeyPoints[i].pt.y = keyPoints[i].y;
        opencvKeyPoints[i].size = keyPoints[i].size;
        opencvKeyPoints[i].angle = keyPoints[i].angle;
    }

    /**draw it*/
    cv::Mat output;
    drawKeypoints(grayImage, opencvKeyPoints, output, Scalar::all(-1), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    imshow("drawKeyPoint", output);

    cvWaitKey(0);
}

/**
 * draw the keypoint on the same pic by SiftySifty and OpenCV
 * @param path
 */
void drawKeyPointCmpToOpenCV(string path)
{
    /**read the pic*/
    cv::Mat originImage = imread(path);
    cv::Mat grayImage;
    cvtColor(originImage, grayImage, CV_RGB2GRAY);

    vector<SiftySifty::KeyPoint> keyPoints;
    SiftySifty::sift(grayImage.data, grayImage.cols, grayImage.rows, keyPoints);

    vector<cv::KeyPoint> opencvKeyPoints1(keyPoints.size());
    for (int i = 0; i < keyPoints.size(); ++i)
    {
        opencvKeyPoints1[i].pt.x  = keyPoints[i].x;
        opencvKeyPoints1[i].pt.y  = keyPoints[i].y;
        opencvKeyPoints1[i].size  = keyPoints[i].size;
        opencvKeyPoints1[i].angle = keyPoints[i].angle;
    }

    cv::Mat output1;
    drawKeypoints(grayImage, opencvKeyPoints1, output1, Scalar::all(-1), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    imshow("SiftySifty", output1);

    vector<cv::KeyPoint> opencvKeyPoints2;
    Ptr<Feature2D> f2d = xfeatures2d::SIFT::create();
    f2d->detect(grayImage, opencvKeyPoints2);

    cv::Mat output2;
    drawKeypoints(grayImage, opencvKeyPoints2, output2, Scalar::all(-1), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    imshow("OpenCV", output2);

    cvWaitKey(0);
}

/**
 * math the keypoint that extracted by SiftySifty
 * @param path
 */
void matchKeyPoint(string path1, string path2)
{
    /**read the pic*/
    cv::Mat originImage1 = imread(path1);
    cv::Mat grayImage1;
    cvtColor(originImage1, grayImage1, CV_RGB2GRAY);

    vector<SiftySifty::KeyPoint> keyPoints1;
    SiftySifty::sift(grayImage1.data, grayImage1.cols, grayImage1.rows, keyPoints1);

    vector<cv::KeyPoint> opencvKeyPoints1(keyPoints1.size());
    cv::Mat ds1(keyPoints1.size(), 128, CV_32F);
    for (int i = 0; i < keyPoints1.size(); ++i)
    {
        opencvKeyPoints1[i].pt.x = keyPoints1[i].x;
        opencvKeyPoints1[i].pt.y = keyPoints1[i].y;
        opencvKeyPoints1[i].size = keyPoints1[i].size;
        opencvKeyPoints1[i].angle = keyPoints1[i].angle;

        memcpy((float *) ds1.data + i * 128, keyPoints1[i].descriptor, sizeof(float) * 128);
    }

    /**read the pic*/
    cv::Mat originImage2 = imread(path2);
    cv::Mat grayImage2;
    cvtColor(originImage2, grayImage2, CV_RGB2GRAY);

    vector<SiftySifty::KeyPoint> keyPoints2;
    SiftySifty::sift(grayImage2.data, grayImage2.cols, grayImage2.rows, keyPoints2);

    vector<cv::KeyPoint> opencvKeyPoints2(keyPoints2.size());
    cv::Mat ds2(keyPoints2.size(), 128, CV_32F);
    for (int i = 0; i < keyPoints2.size(); ++i)
    {
        opencvKeyPoints2[i].pt.x = keyPoints2[i].x;
        opencvKeyPoints2[i].pt.y = keyPoints2[i].y;
        opencvKeyPoints2[i].size = keyPoints2[i].size;
        opencvKeyPoints2[i].angle = keyPoints2[i].angle;

        memcpy((float *) ds2.data + i * 128, keyPoints2[i].descriptor, sizeof(float) * 128);
    }

    BFMatcher matcher;
    vector<DMatch> matches;
    matcher.match(ds1, ds2, matches);

    cv::Mat img_matches;
    drawMatches(grayImage1, opencvKeyPoints1, grayImage2, opencvKeyPoints2, matches, img_matches);
    imshow("SiftySifty match SiftySifty", img_matches);

    cvWaitKey(0);
}

/**
 * match the siftysifty keypoint with opencv keypoint
 * @param path
 */
void matchKeyPointSiftySiftyWithOpenCV(string path1, string path2)
{
    /**read the pic*/
    cv::Mat originImage1 = imread(path1);
    cv::Mat grayImage1;
    cvtColor(originImage1, grayImage1, CV_RGB2GRAY);

    vector<SiftySifty::KeyPoint> keyPoints1;
    SiftySifty::sift(grayImage1.data, grayImage1.cols, grayImage1.rows, keyPoints1);

    vector<cv::KeyPoint> opencvKeyPoints1(keyPoints1.size());
    cv::Mat ds1(keyPoints1.size(), 128, CV_32F);
    for (int i = 0; i < keyPoints1.size(); ++i)
    {
        opencvKeyPoints1[i].pt.x = keyPoints1[i].x;
        opencvKeyPoints1[i].pt.y = keyPoints1[i].y;
        opencvKeyPoints1[i].size = keyPoints1[i].size;
        opencvKeyPoints1[i].angle = keyPoints1[i].angle;

        memcpy((float *) ds1.data + i * 128, keyPoints1[i].descriptor, sizeof(float) * 128);
    }

    /**read the pic*/
    cv::Mat originImage2 = imread(path2);
    cv::Mat grayImage2;
    cvtColor(originImage2, grayImage2, CV_RGB2GRAY);
    vector<cv::KeyPoint> opencvKeyPoints2;
    Ptr<Feature2D> f2d = xfeatures2d::SIFT::create();
    cv::Mat ds2;
    f2d->detectAndCompute(grayImage2, noArray(), opencvKeyPoints2, ds2);

    BFMatcher matcher;
    vector<DMatch> matches;
    matcher.match(ds1, ds2, matches);

    cv::Mat img_matches;
    drawMatches(grayImage1, opencvKeyPoints1, grayImage2, opencvKeyPoints2, matches, img_matches);
    imshow("SiftySifty match OpenCV", img_matches);

    cvWaitKey(0);
}

/**
 * test the speed of SiftySifty and OpenCV
 * @param path
 */
void testSpeedSiftySiftyAndOpenCV(string path) {
    /**read the pic*/
    cv::Mat originImage = imread(path);
    cv::Mat grayImage;
    cvtColor(originImage, grayImage, CV_RGB2GRAY);

    long total = 0;
    for (int i = 0; i < 100; ++i) {
        long t1 = getCurrentTime();

        vector<SiftySifty::KeyPoint> keyPoints;
        sift(grayImage.data, grayImage.cols, grayImage.rows, keyPoints);

        long t2 = getCurrentTime();

        long cur = t2 - t1;
        total += cur;

        cout << "SiftySifty, time:" << (i + 1) << ", cost:" << cur << "ms" << endl;
    }

    float siftysiftyCost = 1.0 * total / 100;

    total = 0;
    for (int i = 0; i < 100; ++i) {
        long t1 = getCurrentTime();

        vector<cv::KeyPoint> opencvKeyPoints;
        Ptr<Feature2D> f2d = xfeatures2d::SIFT::create();
        cv::Mat ds;
        f2d->detectAndCompute(grayImage, noArray(), opencvKeyPoints, ds);

        long t2 = getCurrentTime();

        long cur = t2 - t1;
        total += cur;

        cout << "OpenCV, time:" << (i + 1) << ", cost:" << cur << "ms" << endl;
    }

    float opencvCost = 1.0 * total / 100;

    cout << "SiftySifty cost time(the average of 100 times):" << siftysiftyCost << "ms" << endl;
    cout << "OpenCV cost time(the average of 100 times):" << opencvCost << "ms" << endl;
}

}

















