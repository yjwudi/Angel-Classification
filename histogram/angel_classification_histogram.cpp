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

FILE *fp, *fp1;
Mat ground_img;
MatND ground_hist;
Mat open_img, show_img;
int height, width;
bool points[3500][2500];
vector< pair<int, int> > mid_vec_up[10], mid_vec_down[10];
vector<int> mid_y_up[10], mid_y_down[10];
void elimination();
int check_left(int x, int y);
int check_right(int x, int y);
int calc_zero(int x, int y);
void Cut_img(Mat src_img);
bool get_distance(Mat img);
void insert_midpoint_up(int a, int b, int y);
void insert_midpoint_down(int a, int b, int y);
void get_classification(vector<int> ans_up, vector<int> ans_down);
int main()
{
    /*
    get the points between gum and teeth
    */
    Mat img = imread("/home/yjwudi/cv_code/pic/test/test.jpg");
    cvtColor(img, img, CV_BGR2GRAY);
    vector<int> compression_params;
    compression_params.push_back(IMWRITE_JPEG_QUALITY);
    compression_params.push_back(9);

    fp = fopen("tu_gray_points", "w");
    vec.clear();

    Cut_img(img);
    cout << vec.size() << endl;
    for(int i = 0; i < vec.size(); i++)
    {
        circle(img, vec[i], 1, Scalar( 0, 0, 255 ), -1, 8);
    }
    imwrite("tu_gray_sum.jpg", img, compression_params);
    fclose(fp);
    
    /*
    process the points and eliminate the noise
    */
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
    elimination();


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
    imwrite("tu_eli.jpg", show_img, compression_params);
    fclose(fp);
    fclose(fp1);

    /*
    judge
    */
    FILE *out_f = fopen("output", "w");
    pair<int, int> pii;
    open_img = imread("/home/yjwudi/tpic/tooth/third_gray.jpg");
    fp = fopen("/home/yjwudi/cv_code/points/eli_third__points", "r");
    if(fp == NULL)
    {
        cout << "hehe\n";
        return 0;
    }
    height = open_img.rows;
    width = open_img.cols;
    memset(points, 0, sizeof(points));
    int x, y, z;
    int lower, upper, left_, right_;
    lower = height/4, upper = 3*lower;
    left_ = 250, right_ = 0.8*width;

    while(fscanf(fp, "%d %d", &x, &y)!=EOF)
    {
        if((x>=left_ && x<=right_) && (y>=lower&&y<=upper))
            points[x][y] = true;
    }

    for(i = 0; i < 10; i++)
    {
        mid_y_up[i].clear();
        mid_vec_up[i].clear();
        mid_y_up[i].push_back(-1);
        mid_vec_up[i].push_back(make_pair(0,500*i*2));
        mid_y_down[i].clear();
        mid_vec_down[i].clear();
        mid_y_down[i].push_back(-1);
        mid_vec_down[i].push_back(make_pair(0,500*i*2));
    }

    int left, right;
    for(y = 2*lower; y > lower; y--)
    {
        for(x = left_; x < right_; x++)
        {
            if(points[x][y])
            {
                left = check_left(x, y);
                right = check_right(x, y);
                if(right <= left+2 || left >= 10)
                    continue;
                z = x;
                while(points[z][y])
                    z++;
                x = z;

                while(1)
                {
                    if(z > right_)
                        break;
                    if(points[z][y])
                    {
                        left = check_left(z, y);
                        right = check_right(z, y);
                        if(left > right+2)
                            break;
                    }
                    z++;
                }
                if(z>right_ || z-x>800)
                    break;
                insert_midpoint_up(x, z, y);
                x = z;
                while(points[x][y])
                    x++;
            }
        }
    }
    for(y = 2.5*lower; y < upper; y++)
    {
        for(x = left_; x < right_; x++)
        {
            if(points[x][y])
            {
                left = check_left(x, y);
                right = check_right(x, y);
                if(right <= left+2 || left >= 10)
                    continue;
                z = x;
                while(points[z][y])
                    z++;
                x = z;

                while(1)
                {
                    if(z > right_)
                        break;
                    if(points[z][y])
                    {
                        left = check_left(z, y);
                        right = check_right(z, y);
                        if(left > right+2)
                            break;
                    }
                    z++;
                }
                if(z>right_)
                    break;
                if(z-x>400)
                    continue;
                insert_midpoint_down(x, z, y);
                x = z;
                while(points[x][y])
                    x++;
            }
        }
    }
    vector<int> ans_up, ans_down;
    for(i = 1; i < 6; i++)
    {
        int mid_sum = 0;
        for(j = 1; j < mid_vec_up[i].size(); j++)
        {
            pii = mid_vec_up[i][j];
            mid_sum += (pii.first+pii.second)/2;
            circle(open_img, Point(pii.first, mid_y_up[i][j]), 1, Scalar( 0, 255, 0 ), 5, 8);
            circle(open_img, Point(pii.second, mid_y_up[i][j]), 1, Scalar( 255, 0, 0 ), 5, 8);
            fprintf(out_f,"(%d %d %d %d) ", pii.first, pii.second, (pii.first+pii.second)/2, mid_y_up[i][j]);
        }
        ans_up.push_back(mid_sum/mid_vec_up[i].size());
        fprintf(out_f, "%d\n", mid_sum/mid_vec_up[i].size());
    }
    fprintf(out_f, "\n\ndown:\n");
    for(i = 1; i < 6; i++)
    {
        int mid_sum = 0;
        for(j = 1; j < mid_vec_down[i].size(); j++)
        {
            pii = mid_vec_down[i][j];
            mid_sum += (pii.first+pii.second)/2;
            circle(open_img, Point(pii.first, mid_y_down[i][j]), 1, Scalar( 0, 255, 0 ), 5, 8);
            circle(open_img, Point(pii.second, mid_y_down[i][j]), 1, Scalar( 255, 0, 0 ), 5, 8);
            fprintf(out_f,"(%d %d %d %d) ", pii.first, pii.second, (pii.first+pii.second)/2, mid_y_down[i][j]);
        }
        ans_down.push_back(mid_sum/mid_vec_down[i].size());
        fprintf(out_f, "%d\n", mid_sum/mid_vec_down[i].size());
    }
    cout << "ans up:\n";
    for(i = 0; i < ans_up.size(); i++)
        cout << ans_up[i] << " ";
    cout << endl;
    cout << "ans down:\n";
    for(i = 0; i < ans_down.size(); i++)
        cout << ans_down[i] << " ";
    cout << endl;
    get_classification(ans_up, ans_down);


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
void elimination()
{
    int x, y, left, right;

    for(x=0; x<3500; x++)
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
void get_classification(vector<int> ans_up, vector<int> ans_down)
{
    vector<int> up_vec, down_vec;
    up_vec.clear();
    down_vec.clear();
    int i, j;
    for(i = 0; i < ans_up.size(); i++)
        if(ans_up[i])
        up_vec.push_back(ans_up[i]);
    for(i =0; i < ans_down.size(); i++)
        if(ans_down[i])
        down_vec.push_back(ans_down[i]);
    if(up_vec.size()>=3 && down_vec.size()>=3)
    {
        int up_point, down_point;
        int general_len, molar_len;
        int up_len1=up_vec[1]-up_vec[0], up_len2=up_vec[2]-up_vec[1];
        int down_len1=down_vec[1]-down_vec[0], down_len2=down_vec[2]-down_vec[1];
        if(up_len1 >= 600 || up_len2 >= 600)
        {
            if(up_len2 >= 600)
            {
                up_point = up_vec[2];
                general_len = 2*up_len2/3;
            }
            else
            {
                up_point = up_vec[1];
                general_len = up_len2;
                if(up_point > 1500)
                {
                    up_point = up_vec[0];
                    general_len = up_vec[2]-up_vec[1];
                }
            }
        }
        else
        {
            if(up_len1 > 3*up_len2/2 || abs(up_len1-up_len2)>200)
            {
                up_point = up_vec[1];
                general_len = up_len2;
            }
            else
            {
                up_point = up_vec[0];
                general_len = (up_len1+up_len2)/2;
            }

        }
        if(down_len1 >= 600 || down_len2 >= 600)
        {
            if(down_len2 >= 600)
            {
                down_point = down_vec[2];
            }
            else
            {
                down_point = down_vec[1];
                if(down_point > 1500)
                {
                    down_point = down_vec[0];
                }
            }
        }
        else
        {
            if(down_len1 > 2*down_len2/3 || abs(down_len1-down_len2)>200)
            {
                down_point = down_vec[1];
            }
            else
            {
                down_point = down_vec[0];
            }

        }
        molar_len = 3*general_len/2;
        Debug(down_point);
        Debug(up_point);
        Debug(molar_len);
        if(down_point-up_point>min(molar_len/10,100) && down_point-up_point<=3*molar_len/4)
        {
            cout << "安氏一类\n";
            return ;
        }
        else if(up_point >= down_point || down_point-up_point<=min(molar_len/10,100))
        {
            cout << "安氏二类\n";
            return ;
        }
        else
        {
            cout << "安氏三类\n";
            return ;
        }
    }
    srand((unsigned) time(NULL));
    int number = rand()%2;
    if(number)
    {
        cout << "安氏二类\n";
    }
    else
    {
        cout << "安氏一类\n";
    }
}
void insert_midpoint_up(int a, int b, int y)
{
    int i, j;
    int m = (a+b)/2;
    pair<int, int> pii;
    int dis, min_dis = 100000, min_pos;
    for(i = 0; i < 10; i++)
    {
        dis = 0;
        for(j = 0; j < mid_vec_up[i].size(); j++)
        {
            pii = mid_vec_up[i][j];
            dis += abs((pii.first+pii.second)/2-m);
        }
        if(dis/mid_vec_up[i].size() < min_dis)
        {
            min_dis = dis/mid_vec_up[i].size();
            min_pos = i;
        }
    }
    mid_vec_up[min_pos].push_back(make_pair(a, b));
    mid_y_up[min_pos].push_back(y);
}
void insert_midpoint_down(int a, int b, int y)
{
    int i, j;
    int m = (a+b)/2;
    pair<int, int> pii;
    int dis, min_dis = 100000, min_pos;
    for(i = 0; i < 10; i++)
    {
        dis = 0;
        for(j = 0; j < mid_vec_down[i].size(); j++)
        {
            pii = mid_vec_down[i][j];
            dis += abs((pii.first+pii.second)/2-m);
        }
        if(dis/mid_vec_down[i].size() < min_dis)
        {
            min_dis = dis/mid_vec_down[i].size();
            min_pos = i;
        }
    }
    mid_vec_down[min_pos].push_back(make_pair(a, b));
    mid_y_down[min_pos].push_back(y);
}

