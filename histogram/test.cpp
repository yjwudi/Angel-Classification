#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <vector>
#include <cstring>
using namespace std;
using namespace cv;



int kase = 0;
bool points[3500][2500];
Mat img = imread("/home/yjwudi/face++/face/zjy/test/test.jpg");
//Mat img = imread("/home/yjwudi/tpic/tooth/third_gray.jpg");
FILE *fp1, *fp2;
void get_hist(Mat roi_img);
static void onMouse( int event, int x, int y, int, void* );
void rgb2hsv(int r, int g, int b);
int main()
{
/*
    FILE *fp = fopen("/home/yjwudi/tpic/tooth/points/eli_test__points", "r");

    memset(points, 0, sizeof(points));
    int x, y, z;
    while(fscanf(fp, "%d %d", &x, &y)!=EOF)
    {
        if(x>=250 && x<=1000 && ((y>=690&&y<=1200)||(y>=1350&&y<=1650)))
            points[x][y] = true;
        if(x>=1001 && x<=2800 && ((y>=550&&y<=1200)||(y>=1350&&y<=1880)))
            points[x][y] = true;
    }

    for(x = 0; x < 3500; x++)
    {
        y = 1050;
        //for(y = 0; y < 2500; y++)
        {
            if(points[x][y])
                circle(img, Point(x, y), 1, Scalar( 0, 0, 255 ), 0.5, 8);
        }
    }

*/
    Mat img = imread("hehe.jpg");
    Mat gray_img;
    cvtColor(img, gray_img, CV_BGR2GRAY);
    vector<int> compression_params;
    compression_params.push_back(IMWRITE_JPEG_QUALITY);
    compression_params.push_back(9);
    imwrite("hehe_gray.jpg", gray_img, compression_params);


    namedWindow("src_img", WINDOW_NORMAL);
    imshow("src_img",img);


    cout << img.rows << " " << img.cols << endl;
    setMouseCallback("src_img", onMouse, 0);

    //fclose(fp1);
    //fclose(fp);
    waitKey();
    return 0;
}
void get_hist(Mat roi_img)
{
    int channels = 0;
    //然后是配置输出的结果存储的空间 ，用MatND类型来存储结果
    MatND dstHist;
    //接下来是直方图的每一个维度的 柱条的数目（就是将数值分组，共有多少组）
    int histSize[] = { 32 };
    float midRanges[] = { 0, 256 };
    const float *ranges[] = { midRanges };

    calcHist(&roi_img, 1, &channels, Mat(), dstHist, 1, histSize, ranges, true, false);

    //calcHist  函数调用结束后，dstHist变量中将储存了 直方图的信息  用dstHist的模版函数 at<Type>(i)得到第i个柱条的值
    //at<Type>(i, j)得到第i个并且第j个柱条的值

    //因为任何一个图像的某个像素的总个数，都有可能会有很多，会超出所定义的图像的尺寸，针对这种情况，先对个数进行范围的限制
    //先用 minMaxLoc函数来得到计算直方图后的像素的最大个数
    double g_dHistMaxValue;
    minMaxLoc(dstHist, 0, &g_dHistMaxValue, 0, 0);
    //将像素的个数整合到 图像的最大范围内
    //遍历直方图得到的数据
    RNG &rng = theRNG();
    int sum = 0;
    for (int i = 0; i < 32; i++)
    {
        int value = cvRound(dstHist.at<float>(i) * 256 * 0.9 / g_dHistMaxValue);
        cout << value << " ";
        sum += value;
        //
    }
    cout << endl << sum << endl;
}
static void onMouse( int event, int x, int y, int, void* )
{
    //cout << x << " " << y << endl;
    switch(event)
    {
    case CV_EVENT_LBUTTONUP:

        int r, g, b;
        r = (int)img.at<Vec3b>(x,y)[2];
        g = (int)img.at<Vec3b>(x,y)[1];
        b = (int)img.at<Vec3b>(x,y)[0];
        //rgb2hsv(r,g,b);
        /*
        cout << x << " " << y << endl;
        cout << (int)img.at<Vec3b>(x,y)[0] << " ";
        cout << (int)img.at<Vec3b>(x,y)[1] << " ";
        cout << (int)img.at<Vec3b>(x,y)[2] << endl;
        int sum = 0;
        sum = (int)img.at<Vec3b>(x,y)[0]+(int)img.at<Vec3b>(x,y)[1]+(int)img.at<Vec3b>(x,y)[2];
        cout << "sum: " << sum << endl;
        */
        fp1 = fopen("/home/yjwudi/tpic/tooth/white_rgb", "a");
        int sum1 = 0, sum2 = 0, sum3 = 0;
        int i, j;
        for(i = x; i < x+20; i++)
        {
            for(j = y; j < y+20; j++)
            {
                sum1 += (int)img.at<Vec3b>(i,j)[0];
                sum2 += (int)img.at<Vec3b>(i,j)[1];
                sum3 += (int)img.at<Vec3b>(i,j)[2];
            }
        }

        //cout << sum1 << " " << sum2 << " " << sum3 << endl;
        fprintf(fp1, "%d %d %d %d\n", sum3, sum2, sum1, sum1+sum2+sum3);
        fclose(fp1);

        Mat roi_img;
        cout << x << " " << y << endl;
        Rect rect(x,y,30,30);
        img(rect).copyTo(roi_img);
        get_hist(roi_img);
    }
}
void rgb2hsv(int r, int g, int b)
{
//cout << r << " " << g << " " << b << endl;
    double h, s, v;
    v = max(r, max(g, b));
    if(v == 0)
        s = 0;
    else
        s = (v-min(r, min(g,b)))/v;
    if((int)v == r)
    {
        h = 60*(g-b)/(v-min(r, min(g,b)));
    }
    else if((int)v == g)
    {
        h = 120+60*(b-r)/(v-min(r, min(g,b)));
    }
    else if((int)v == b)
    {
        h = 240+60*(r-g)/(v-min(r, min(g,b)));
    }
    if(h < 0)
        h += 360;
    cout << h << " " << s << " " << v << endl;
}











