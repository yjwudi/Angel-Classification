#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <cstdio>
#include <vector>
using namespace std;
using namespace cv;

vector<Point> vec;
//剪切图片为m * n 块

FILE *fp;
Mat ground_img;
MatND ground_hist;
void Cut_img(Mat src_img);
bool get_distance(Mat img);
int main()
{
    Mat img;// = imread("/home/yjwudi/cv_code/pic/test/test.jpg");
    //Mat gray_img;
    //cvtColor(img, gray_img, CV_BGR2GRAY);
    vector<int> compression_params;
    compression_params.push_back(IMWRITE_JPEG_QUALITY);
    compression_params.push_back(9);
    //imwrite("/home/yjwudi/cv_code/pic/test/test_gray.jpg", gray_img, compression_params);

    img = imread("tu_gray.jpg");

    fp = fopen("tu_gray_points", "w");
    vec.clear();

    Cut_img(img);
    cout << vec.size() << endl;
    for(int i = 0; i < vec.size(); i++)
    {
        circle(img, vec[i], 1, Scalar( 0, 0, 255 ), -1, 8);
    }
    imwrite("tu_gray_sum.jpg", img, compression_params);
    namedWindow("src img", WINDOW_NORMAL);
    imshow("src img",img);
    fclose(fp);

    waitKey();
    return 0;

}
bool get_distance(Mat img)
{
    int channels = 0;
    MatND dstHist;
    int histSize[] = { 32 };
    float midRanges[] = { 0, 256 };
    const float *ranges[] = { midRanges };

    calcHist(&img, 1, &channels, Mat(), dstHist, 1, histSize, ranges, true, false);

    double g_dHistMaxValue;
    minMaxLoc(dstHist, 0, &g_dHistMaxValue, 0, 0);
    //将像素的个数整合到 图像的最大范围内
    //遍历直方图得到的数据
    RNG &rng = theRNG();
    int sum = 0;
    for (int i = 0; i < 32; i++)
    {
        int value = cvRound(dstHist.at<float>(i) * 256 * 0.9 / g_dHistMaxValue);
        sum += value;
    }
    if(sum > 600)
    {
        return true;
    }
    return false;
}
void Cut_img(Mat src_img)
{
    int height = src_img.rows;
    int width  = src_img.cols;
    cout << "rows: " << src_img.rows << endl;
    cout << "cols: " << src_img.cols << endl;

    int ceil_height = 50;
    int ceil_width  = 50;

    Mat roi_img,tmp_img;

    for(int i = height/4; i < height-height/4; i++)
        for(int j = 250; j < 0.8*width; j++)
        {
            if(i < height/4+200 && j < 1000)
                continue;
            //x(j)是横向的坐标，y(i)是纵向的坐标
            Rect rect(j,i,ceil_width,ceil_height);
            //fprintf(fp, "%d %d ", i, j);
            //cout << i << " " << j << " ";
            src_img(rect).copyTo(roi_img);
            if(get_distance(roi_img))
            {
                Point p(j+25, i+25);
                fprintf(fp, "%d %d\n", j+25, i+25);
                vec.push_back(p);
            }
        }
}



