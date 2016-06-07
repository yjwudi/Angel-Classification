#include "slic.h"
#include "slic.cpp"
#include <cstdlib>
#include <cstdio>

using namespace std;

SLIC slic;
static void onMouse( int event, int x, int y, int, void* );
int main(int argc, char** argv)
{
	if (argc != 3) {
		printf("usage: test_slic <filename> <number of superpixels>\n");
		exit(-1);
	}

	cv::Mat img, result;
	
	img = imread(argv[1]);
	pyrDown( img, img, Size( img.cols/2, img.rows/2 ) );
	int numSuperpixel = atoi(argv[2]);

	
	slic.GenerateSuperpixels(img, numSuperpixel);
	return slic.GetImgWithContours(cv::Scalar(0, 0, 255));
	/*
	if (img.channels() == 3) 
		result = slic.GetImgWithContours(cv::Scalar(0, 0, 255));
	else
		result = slic.GetImgWithContours(cv::Scalar(128));

	namedWindow("src_img", WINDOW_NORMAL);
    imshow("src_img",result);
    setMouseCallback("src_img", onMouse, 0);
    waitKey();
	*/
	/*
	slic.TestEdge(img.rows, img.cols);
	namedWindow("src_img", WINDOW_NORMAL);
    imshow("src_img",img);
    setMouseCallback("src_img", onMouse, 0);
    waitKey();
    */
    
	
	//cv::imwrite("result_1000_green.jpg", result);
}

static void onMouse( int event, int x, int y, int, void* )
{
    switch(event)
    {
        case EVENT_LBUTTONUP:
            cout << "lbuttonup: " << x << " " << y << endl;
            ///slic.GetRegionInfo(x, y);
    }
}
