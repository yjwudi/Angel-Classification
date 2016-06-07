#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <cstdio>
#include <vector>
#include <cstring>
using namespace std;
using namespace cv;

vector<Point> vec;
//剪切图片为m * n 块

FILE *fp, *fp1;
Mat open_img, show_img;
int height, width;
bool points[3500][2500];
void Cut_img(Mat src_img);
bool get_distance(Mat img);
void elimination();
int check_left(int x, int y);
int check_right(int x, int y);
int calc_zero(int x, int y);
static void onMouse( int event, int x, int y, int, void* );
int main()
{
    int i, j;
    open_img = imread("tu_gray.jpg");
    show_img = imread("tu_gray.jpg");
    fp = fopen("tu_gray_points", "r");
    fp1 = fopen("eli_tu__points", "w");
    height = open_img.rows;
    width = open_img.cols;
    int lower, upper, left_, right_;
    lower = height/4, upper = 3*lower;
    left_ = 250, right_ = 0.8*width;
    memset(points, 0, sizeof(points));
    int x, y;
    while(fscanf(fp, "%d %d", &x, &y)!=EOF)
    {
        if((x>=left_ && x<=right_) && (y>=lower&&y<=upper))
                {
                    points[x][y] = true;
                    circle(open_img, Point(x, y), 1, Scalar( 0, 0, 255 ), -1, 8);
                }

    }
    //elimination();
    namedWindow("src_img", WINDOW_NORMAL);
    imshow("src_img",open_img);
    setMouseCallback("src_img", onMouse, 0);
    waitKey();


    for(i = 0; i < 3500; i++)
    {
        for(j = 0; j < 2500; j++)
        {
            if(points[i][j])
            {
                fprintf(fp1, "%d %d\n", i, j);
                circle(show_img, Point(i, j), 1, Scalar( 0, 0, 255 ), -1, 8);
            }
        }
    }
    namedWindow("src", WINDOW_NORMAL);
    imshow("src",show_img);

    vector<int> compression_params;
    compression_params.push_back(IMWRITE_JPEG_QUALITY);
    compression_params.push_back(9);
    imwrite("tu_eli.jpg", show_img, compression_params);


    waitKey();
    fclose(fp);
    fclose(fp1);
    return 0;
}
static void onMouse( int event, int x, int y, int, void* )
{
    
    switch(event)
    {
    case CV_EVENT_LBUTTONUP:
    cout << x << " " << y << endl;
    int i, j;
    for(i = x; i < x+100; i++)
    {
        for(j = y; j < y+100; j++)
        {
            points[i][j] = false;
        }
    }
    }
}
void elimination()
{
    int x, y, left, right;

    for(x=0; x<3500; x++)//y
    {
        for(y=0; y<2500; y++)
        {
            if(points[x][y])
            {
                //浅色是1，深色是2
                left = check_left(x, y);
                right = check_right(x, y);
                if(abs(left - right) < 2)
                {
                    points[x][y] = false;
                }
            }
        }
    }
}
int check_left(int x, int y)
{
    int i, j, k, sum;
    int zero_num1 = 0, zero_num2 = 0;
    for(i=x-10; i>50; i--)
    {
        if(!points[i][y])
        {
            i = i-30;
            if(i > 10)
            {
                zero_num1 = calc_zero(i, y);
            }
            i = i-30;
            if(i > 10)
            {
                zero_num2 = calc_zero(i, y);
            }
            break;
        }
    }
    return (zero_num1+zero_num2)/2;
}
int check_right(int x, int y)
{
    int i, j, k, sum;
    int zero_num1 = 0, zero_num2 = 0;
    for(i=x+10; i<0.8*width; i++)
    {
        if(!points[i][y])
        {
            if(i<0.8*width)
                zero_num1 = calc_zero(i, y);
            i = i+30;
            if(i<0.8*width)
                zero_num2 = calc_zero(i, y);
            break;
        }
    }
    return (zero_num1+zero_num2)/2;
}

int calc_zero(int x, int y)
{
    Mat img;
    Rect rect(x,y,30,30);
    open_img(rect).copyTo(img);
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
    for (int i = 31; i >= 0; i--)
    {
        int value = cvRound(dstHist.at<float>(i) * 256 * 0.9 / g_dHistMaxValue);
        if(value == 0)
            sum++;
        else
            break;
    }
    return sum;
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
    if(sum > 750)
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

    for(int i = 0; i < height-50; i++)
        for(int j = 0; j < width-50; j++)
        {
            //x(j)是横向的坐标，y(i)是纵向的坐标
            Rect rect(j,i,ceil_width,ceil_height);
            src_img(rect).copyTo(roi_img);
            if(get_distance(roi_img))
            {
                Point p(j+25, i+25);
                fprintf(fp, "%d %d\n", j+25, i+25);//都是patch的中心点
                vec.push_back(p);
            }
        }
}



