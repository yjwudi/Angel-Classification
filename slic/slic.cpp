// SLIC.cpp: implementation of the SLIC class.
//===========================================================================
// This code implements the zero parameter superpixel segmentation technique
// described in:
//
//
//
// "SLIC Superpixels Compared to State-of-the-art Superpixel Methods"
//
// Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal Fua,
// and Sabine Susstrunk,
//
// IEEE TPAMI, Volume 34, Issue 11, Pages 2274-2282, November 2012.
//
//
//===========================================================================
// Copyright (c) 2013 Radhakrishna Achanta.
//
// For commercial use please contact the author:
//
// Email: firstname.lastname@epfl.ch
//===========================================================================
//
// Modified by nipan
// Email: nipan1988@gmail.com
//

#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include "slic.h"

// For superpixels
const int dx4[4] = {-1,  0,  1,  0};
const int dy4[4] = { 0, -1,  0,  1};
//const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
//const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

// For supervoxels
const int dx10[10] = {-1,  0,  1,  0, -1,  1,  1, -1,  0, 0};
const int dy10[10] = { 0, -1,  0,  1, -1, -1,  1,  1,  0, 0};
const int dz10[10] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, 1};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SLIC::SLIC()
{
    m_lvec = NULL;
    m_avec = NULL;
    m_bvec = NULL;

    m_lvecvec = NULL;
    m_avecvec = NULL;
    m_bvecvec = NULL;

    bufferGray = NULL;
    bufferRGB = NULL;
}

SLIC::~SLIC()
{
    if(m_lvec) delete [] m_lvec;
    if(m_avec) delete [] m_avec;
    if(m_bvec) delete [] m_bvec;


    if(m_lvecvec)
    {
        for( int d = 0; d < m_depth; d++ ) delete [] m_lvecvec[d];
        delete [] m_lvecvec;
    }
    if(m_avecvec)
    {
        for( int d = 0; d < m_depth; d++ ) delete [] m_avecvec[d];
        delete [] m_avecvec;
    }
    if(m_bvecvec)
    {
        for( int d = 0; d < m_depth; d++ ) delete [] m_bvecvec[d];
        delete [] m_bvecvec;
    }

    if (bufferGray)
    {
        delete [] bufferGray;
    }

    if (bufferRGB)
    {
        delete [] bufferRGB;
    }

    if (label)
    {
        delete [] label;
    }
}

//==============================================================================
///	RGB2XYZ
///
/// sRGB (D65 illuninant assumption) to XYZ conversion
//==============================================================================
void SLIC::RGB2XYZ(
    const int&		sR,
    const int&		sG,
    const int&		sB,
    double&			X,
    double&			Y,
    double&			Z)
{
    double R = sR/255.0;
    double G = sG/255.0;
    double B = sB/255.0;

    double r, g, b;

    if(R <= 0.04045)	r = R/12.92;
    else				r = pow((R+0.055)/1.055,2.4);
    if(G <= 0.04045)	g = G/12.92;
    else				g = pow((G+0.055)/1.055,2.4);
    if(B <= 0.04045)	b = B/12.92;
    else				b = pow((B+0.055)/1.055,2.4);

    X = r*0.4124564 + g*0.3575761 + b*0.1804375;
    Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
    Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
}

//===========================================================================
///	RGB2LAB
//===========================================================================
void SLIC::RGB2LAB(const int& sR, const int& sG, const int& sB, double& lval, double& aval, double& bval)
{
    //------------------------
    // sRGB to XYZ conversion
    //------------------------
    double X, Y, Z;
    RGB2XYZ(sR, sG, sB, X, Y, Z);

    //------------------------
    // XYZ to LAB conversion
    //------------------------
    double epsilon = 0.008856;	//actual CIE standard
    double kappa   = 903.3;		//actual CIE standard

    double Xr = 0.950456;	//reference white
    double Yr = 1.0;		//reference white
    double Zr = 1.088754;	//reference white

    double xr = X/Xr;
    double yr = Y/Yr;
    double zr = Z/Zr;

    double fx, fy, fz;
    if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
    else				fx = (kappa*xr + 16.0)/116.0;
    if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
    else				fy = (kappa*yr + 16.0)/116.0;
    if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
    else				fz = (kappa*zr + 16.0)/116.0;

    lval = 116.0*fy-16.0;
    aval = 500.0*(fx-fy);
    bval = 200.0*(fy-fz);
}

//===========================================================================
///	DoRGBtoLABConversion
///
///	For whole image: overlaoded floating point version
//===========================================================================
void SLIC::DoRGBtoLABConversion(
    const unsigned int*&		ubuff,
    double*&					lvec,
    double*&					avec,
    double*&					bvec)
{
    int sz = m_width*m_height;
    lvec = new double[sz];
    avec = new double[sz];
    bvec = new double[sz];

    for( int j = 0; j < sz; j++ )
    {
        int r = (ubuff[j] >> 16) & 0xFF;
        int g = (ubuff[j] >>  8) & 0xFF;
        int b = (ubuff[j]      ) & 0xFF;

        RGB2LAB( r, g, b, lvec[j], avec[j], bvec[j] );
    }
}

//===========================================================================
///	DoRGBtoLABConversion
///
/// For whole volume
//===========================================================================
void SLIC::DoRGBtoLABConversion(
    const unsigned int**&		ubuff,
    double**&					lvec,
    double**&					avec,
    double**&					bvec)
{
    int sz = m_width*m_height;
    for( int d = 0; d < m_depth; d++ )
    {
        for( int j = 0; j < sz; j++ )
        {
            int r = (ubuff[d][j] >> 16) & 0xFF;
            int g = (ubuff[d][j] >>  8) & 0xFF;
            int b = (ubuff[d][j]      ) & 0xFF;

            RGB2LAB( r, g, b, lvec[d][j], avec[d][j], bvec[d][j] );
        }
    }
}

//=================================================================================
/// DrawContoursAroundSegments
///
/// Internal contour drawing option exists. One only needs to comment the if
/// statement inside the loop that looks at neighbourhood.
//=================================================================================
void SLIC::DrawContoursAroundSegments(
    unsigned int*			ubuff,
    const int*				labels,
    const int&				width,
    const int&				height,
    const cv::Scalar&		color )
{
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    int sz = width*height;


    vector<bool> istaken(sz, false);

    int mainindex(0);
    int i, j, k;
    cv::Scalar green(0, 255, 0);
    for( j = 0; j < height; j++ )
    {
        for( k = 0; k < width; k++ )
        {
            int np(0);
            for( i = 0; i < 8; i++ )
            {
                int x = k + dx8[i];
                int y = j + dy8[i];

                if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                {
                    int index = y*width + x;

                    if( false == istaken[index] )//comment this to obtain internal contours
                    {
                        if( labels[mainindex] != labels[index] ) np++;
                    }
                }
            }
            if( np > 1 )//change to 2 or 3 for thinner lines
            {
                ubuff[mainindex] = 0;
                ubuff[mainindex] |= (int)color.val[2] << 16; // r
                ubuff[mainindex] |= (int)color.val[1] << 8; // g
                ubuff[mainindex] |= (int)color.val[0];
                istaken[mainindex] = true;
                edge_flag[mainindex] = 1;
            }
            else
            {
                edge_flag[mainindex] = 0;
            }
            //把牙齿标记为绿色
            if(label_flag[label[mainindex]] == 1)
            {
                ubuff[mainindex] = 0;
                ubuff[mainindex] |= (int)green.val[2] << 16; // r
                ubuff[mainindex] |= (int)green.val[1] << 8; // g
                ubuff[mainindex] |= (int)green.val[0];
                point_flag[mainindex] = 1;
            }
            else
            {
                point_flag[mainindex] = 0;
            }
            mainindex++;
        }
    }
    /*
    FILE *fp = fopen("dat_edge.txt", "w");
    for(i = 0; i < mainindex; i++)
    {
        fprintf(fp, "%d\n", edge_flag[i]);
    }
    fclose(fp);
    fp = fopen("dat_point.txt", "w");
    for(i = 0; i < mainindex; i++)
    {
        fprintf(fp, "%d\n", point_flag[i]);
    }
    fclose(fp);
    */
    //TestEdge(m_height, m_width);
}
int SLIC::TestEdge(int height, int width)
{
    int i, j, k;
    int sz = height * width;
    //Mat src = imread("/home/yjwudi/Desktop/slic/san.jpg");
    /*
    point_flag = new int[sz];
    edge_flag = new int[sz];
    FILE *fp = fopen("/home/yjwudi/Desktop/slic/dat_edge.txt", "r");

    i = 0;
    while(fscanf(fp, "%d", &edge_flag[i])!=EOF)
        i++;
    fclose(fp);
    fp = fopen("dat_point.txt", "r");
    i = 0;
    while(fscanf(fp, "%d", &point_flag[i])!=EOF)
        i++;
    fclose(fp);
    */
    vector<int> ans_up, ans_down;
    ans_up.clear();
    ans_down.clear();
    int rsum = 50;
    int lower, upper, left_, right_;
    GetLowerUpper(height, width, lower, upper);
    //cout << "lower upper: " << lower << " " << upper << endl;
    //lower = height/4, upper = 3*lower;
    left_ = 250, right_ = 0.8*width;
    vector<double> upmid_vec;
    upmid_vec.clear();
    for(i = 0; i < 20; i++)
    {
        mid_vec_up[i].clear();
        mid_vec_down[i].clear();
    }
    {
        //maxillar
        for(i = (lower+upper)/2; i > lower; i--)
        {
            for(j = left_; j < right_; )
            {
                int index = i*width+j;
                if(edge_flag[index] && point_flag[index-5] == 1 && point_flag[index+5] == 0)
                {
                    int left_pos  = j;
                    int right_pos = 0;
                    for(j+=10; j < right_; j++)
                    {
                        int index2 =  i*width+j;
                        if(edge_flag[index2] && point_flag[index2-5] == 0 && point_flag[index2+5] == 1)
                        {
                            right_pos = j;
                            break;
                        }
                    }
                    if(right_pos != 0 && (right_pos-left_pos)<550)
                    {
                        //InsertMidpoint_Up(left_pos, right_pos, i);
                        //circle(src, Point(left_pos, i), 1, Scalar( 0, 255, 0 ), 5, 8);
                        //circle(src, Point(right_pos, i), 1, Scalar( 255, 0, 0 ), 5, 8);
                        upmid_vec.push_back((left_pos+right_pos)/2);
                    }
                }
                else
                {
                    j++;
                }
            }
        }
        //cluster for maxillar
        sort(upmid_vec.begin(), upmid_vec.end());
        int kind_num = 0;
        int begin_ = 0, move_flag = 1;
        int range_sum = 0, last_mid_value = upmid_vec[0];
        for(i = 0; i < upmid_vec.size(); i++)
        {
            if(move_flag)
            {
                if(upmid_vec[i]-last_mid_value < 20)
                {
                    last_mid_value = upmid_vec[i];
                    range_sum++;
                }
                else
                {
                    if(range_sum >rsum)
                    {
                        for(j = begin_; j < i; j++)
                        {
                            mid_vec_up[kind_num].push_back(upmid_vec[j]);
                        }
                        kind_num++;
                        //cout << kind_num << endl;
                    }
                    move_flag = 0;
                }
            }
            else
            {
                range_sum = 0, last_mid_value = upmid_vec[i];
                begin_ = i, move_flag = 1;
            }
        }
        if(move_flag && range_sum > rsum)
        {
            for(j = begin_; j < i; j++)
            {
                mid_vec_up[kind_num].push_back(upmid_vec[j]);
            }
            kind_num++;
        }
        //delete noise range
        int mid_up_flag[200];
        //average
        int sum_len = 0, ave_len = 0;
        int sum_value = 0;
        //cout << "info:\n";
        for(i = 0; i < kind_num; i++)
        {
            sum_len += mid_vec_up[i].size();
            sum_value = 0;
            for(j = 0; j < mid_vec_up[i].size(); j++)
            {
                sum_value += mid_vec_up[i][j];
            }
            mid_vec_up[i].push_back(sum_value/mid_vec_up[i].size());
            //cout << mid_vec_up[i].size() << " " << mid_vec_up[i][mid_vec_up[i].size()-1] << endl;
        }
        ave_len = sum_len/kind_num;
        for(i = 0; i < kind_num-1; i++)
        {
            int delta = mid_vec_up[i+1][mid_vec_up[i+1].size()-1]-mid_vec_up[i][mid_vec_up[i].size()-1];
            if(mid_vec_up[i].size() < ave_len && delta < 150)
            {
                if(mid_vec_up[i].size()<mid_vec_up[i+1].size())
                {
                    mid_up_flag[i] = 0;
                    mid_up_flag[i+1] = 1;
                }
                else
                {
                    mid_up_flag[i] = 1;
                    mid_up_flag[i+1] = 0;
                }
            }
            else
            {
                mid_up_flag[i] = 1;
                if(i == kind_num-2)
                    mid_up_flag[i+1] = 1;
            }
        }
        //cout << "up flag: ";
        for(i = 0; i < kind_num; i++)
        {
            //cout << mid_up_flag[i] << " ";
            if(mid_up_flag[i])
                ans_up.push_back(mid_vec_up[i][mid_vec_up[i].size()-1]);
        }
        /*cout << "\nans up:\n";
        for(i = 0; i < ans_up.size(); i++)
            cout << ans_up[i] << " ";
        cout << endl;*/
    }
    {
        //mandibular
        vector<double> downmid_vec;
        downmid_vec.clear();
        for(i = (lower+upper)/2; i < upper; i++)
        {
            for(j = left_; j < right_; )
            {
                int index = i*width+j;
                if(edge_flag[index] && point_flag[index-5] == 1 && point_flag[index+5] == 0)
                {
                    int left_pos  = j;
                    int right_pos = 0;
                    for(j+=10; j < right_; j++)
                    {
                        int index2 =  i*width+j;
                        if(edge_flag[index2] && point_flag[index2-5] == 0 && point_flag[index2+5] == 1)
                        {
                            right_pos = j;
                            break;
                        }
                    }
                    if(right_pos != 0 && (right_pos-left_pos)<250)
                    {
                        //InsertMidpoint_Up(left_pos, right_pos, i);
                        //circle(src, Point(left_pos, i), 1, Scalar( 0, 255, 0 ), 5, 8);
                        //circle(src, Point(right_pos, i), 1, Scalar( 255, 0, 0 ), 5, 8);
                        downmid_vec.push_back((left_pos+right_pos)/2);
                    }
                }
                else
                {
                    j++;
                }
            }
        }
        //cluster for maxillar
        sort(downmid_vec.begin(), downmid_vec.end());
        int kind_num = 0;
        int begin_ = 0, move_flag = 1;
        int range_sum = 0, last_mid_value = downmid_vec[0];
        for(i = 0; i < downmid_vec.size(); i++)
        {
            if(move_flag)
            {
                if(downmid_vec[i]-last_mid_value < 20)
                {
                    last_mid_value = downmid_vec[i];
                    range_sum++;
                }
                else
                {
                    if(range_sum > rsum)
                    {
                        for(j = begin_; j < i; j++)
                        {
                            mid_vec_down[kind_num].push_back(downmid_vec[j]);
                        }
                        kind_num++;
                    }
                    move_flag = 0;
                }
            }
            else
            {
                range_sum = 0, last_mid_value = downmid_vec[i];
                begin_ = i, move_flag = 1;
            }
        }
        if(move_flag && range_sum > rsum)
        {
            for(j = begin_; j < i; j++)
            {
                mid_vec_down[kind_num].push_back(downmid_vec[j]);
            }
            kind_num++;
        }

        //delete noise range
        int mid_down_flag[200];
        //average
        int sum_len = 0, ave_len = 0;
        int sum_value = 0;
        //cout << "info:\n";
        for(i = 0; i < kind_num; i++)
        {
            sum_len += mid_vec_down[i].size();
            sum_value = 0;
            for(j = 0; j < mid_vec_down[i].size(); j++)
            {
                sum_value += mid_vec_down[i][j];
            }
            mid_vec_down[i].push_back(sum_value/mid_vec_down[i].size());
            //cout << mid_vec_down[i].size() << " " << mid_vec_down[i][mid_vec_down[i].size()-1] << endl;
        }
        ave_len = sum_len/kind_num;
        for(i = 0; i < kind_num-1; i++)
        {
            int delta = mid_vec_down[i+1][mid_vec_down[i+1].size()-1]-mid_vec_down[i][mid_vec_down[i].size()-1];
            if(mid_vec_down[i].size() < ave_len && delta < 150)
            {
                if(mid_vec_down[i].size()<mid_vec_down[i+1].size())
                {
                    mid_down_flag[i] = 0;
                    mid_down_flag[i+1] = 1;
                }
                else
                {
                    mid_down_flag[i] = 1;
                    mid_down_flag[i+1] = 0;
                }
            }
            else
            {
                mid_down_flag[i] = 1;
                if(i == kind_num-2)
                    mid_down_flag[i+1] = 1;
            }
        }
        FILE *out_f = fopen("downmid.txt", "w");
        for(i = 0; i < downmid_vec.size(); i++)
            fprintf(out_f, "%lf\n", downmid_vec[i]);
        fclose(out_f);
        //cout << "down flag: ";
        for(i = 0; i < kind_num; i++)
        {
            //cout << mid_down_flag[i] << " ";
            if(mid_down_flag[i])
                ans_down.push_back(mid_vec_down[i][mid_vec_down[i].size()-1]);
        }
        /*cout << "\nans down:\n";
        for(i = 0; i < ans_down.size(); i++)
            cout << ans_down[i] << " ";
        cout << endl;*/
    }
    return GetClassification(ans_up, ans_down);

    /*
        FILE *out_f = fopen("output", "w");
        for(i = 0; i < kind_num; i++)
        {
        	for(j = 0; j < mid_vec_up[i].size(); j++)
        		fprintf(out_f,"%d  ", mid_vec_up[i][j]);
        	fprintf(out_f, "\n");
        }
        fclose(out_f);

        FILE *out_f = fopen("output", "w");
        vector<int> ans_up, ans_down;
        for(i = 1; i < 6; i++)
        {
            int mid_sum = 0;
            for(j = 0; j < mid_vec_up[i].size(); j++)
            {
                //pii = mid_vec_up[i][j];
                mid_sum += mid_vec_up[i][j];
                fprintf(out_f,"%d  ", mid_vec_up[i][j]);
                //cout << "(" << pii.first << " " << pii.second << " " << (pii.first+pii.second)/2 << ") ";
            }
            ans_up.push_back(mid_sum/mid_vec_up[i].size());
            fprintf(out_f, "%lf\n", (double)mid_sum/mid_vec_up[i].size());
            //cout << mid_sum/mid_vec[i].size() << endl;
            //fprintf(out_f, "\n\n");
        }
        fprintf(out_f, "\n\ndown:\n");
        for(i = 1; i < 6; i++)
        {
            int mid_sum = 0;
            for(j = 1; j < mid_vec_down[i].size(); j++)
            {
                pii = mid_vec_down[i][j];
                mid_sum += (pii.first+pii.second)/2;
                fprintf(out_f,"(%d %d %d %d) ", pii.first, pii.second, (pii.first+pii.second)/2, mid_y_down[i][j]);
            }
            ans_down.push_back(mid_sum/mid_vec_down[i].size());
            fprintf(out_f, "%d\n", mid_sum/mid_vec_down[i].size());
        }
        //fclose(out_f);

        cout << endl;
        cout << "ans down:\n";
        for(i = 0; i < ans_down.size(); i++)
            cout << ans_down[i] << " ";
        cout << endl;


        FILE *out_f = fopen("upmid.txt", "w");
        for(i = 0; i < upmid_vec.size(); i++)
        	fprintf(out_f, "%lf\n", upmid_vec[i]);
        fclose(out_f);
    
    namedWindow("src_win", WINDOW_NORMAL);
    imshow("src_win", src);
    waitKey();
    */

}
int SLIC::GetClassification(vector<int> ans_up, vector<int> ans_down)
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
            if(down_len1 > 3*down_len2/3 || abs(down_len1-down_len2)>200)
            {
                down_point = down_vec[1];
            }
            else
            {
                down_point = down_vec[0];
            }

        }
        molar_len = 3*general_len/2;
        /*Debug(down_point);
        Debug(up_point);
        Debug(molar_len);*/
        if(down_point-up_point>min(molar_len/10,100) && down_point-up_point<=min(3*molar_len/4, 250))
        {
            //cout << "安氏一类\n";
            return 1;
        }
        else if(up_point >= down_point || down_point-up_point<=min(molar_len/10,100))
        {
            //cout << "安氏二类\n";
            return 2;
        }
        else if(down_point-up_point>=min((int)0.8*molar_len, 350))
        {
            //cout << "安氏三类\n";
            return 3;
        }
    }
    //cout << "can't decide\n";
    return 0;
    /*
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
    */
}
void SLIC::GetLowerUpper(int height, int width, int &lower, int &upper)
{
    int i,j;
    int flagnum, min_flagnum;
    int left_ = 250, right_ = 0.8*width;
    min_flagnum = 1000000;
    for(i = 0; i < height/2; i++)
    {
        flagnum = 0;
        for(j = left_; j < right_; j++)
        {
            int index = i*width+j;
            if(point_flag[index])
            {
                flagnum++;
            }
        }
        if(flagnum < min_flagnum)
        {
            min_flagnum = flagnum;
            lower = i;
        }
    }
    min_flagnum = 1000000;
    for(i = height/2; i < height; i++)
    {
        flagnum = 0;
        for(j = left_; j < right_; j++)
        {
            int index = i*width+j;
            if(point_flag[index])
            {
                flagnum++;
            }
        }
        if(flagnum < min_flagnum)
        {
            min_flagnum = flagnum;
            upper = i;
        }
    }
}
/*
void SLIC::SmartCluster(vector<double> can_vec, vector<double> &ans_vec)
{
    int i, j, k, s, t;
    int min_error_k;
    int len = can_vec.size();
    int center_pos[len+1];
    double dist, min_dist;
    double err_dis, min_err_dis;
    vector<double> ans[10];
    min_err_dis = 10000000;
    for(k = 3; k <= 7; k++)
    {
        ans[k].clear();
        double centers[k];
        err_dis = 0;
        int ite = 0;
        memset(center_pos, -1, sizeof(center_pos));
        for(i = 0; i < k; i++)
        {
            //cout << i*len/(k+10) << endl;
            centers[i] = (can_vec[len-1]-can_vec[0])/(k+2)*i;//can_vec[i*len/(k+10)];
        }
        for(i = 0; i < k; i++)
            cout << centers[i] << "  ";
        cout << endl;
        while(!ite)
        {
            ite = 1;
            for(i = 0; i < len; i++)
            {
                min_dist = 10000000;
                int center_flag = 0;
                for(j = 0; j < k; j++)
                {
                    dist = fabs(can_vec[i]-centers[j]);
                    if(dist < min_dist)
                    {
                        min_dist = dist;
                        center_flag = j;
                    }
                }

                if(i == 1000)
                {
                    cout << "1000: ";
                    cout << can_vec[i] << " " << center_flag << " " << centers[center_flag] << endl;
                }
                if(center_flag!=center_pos[i])
                    ite = 0;
                center_pos[i] = center_flag;
            }
            for(j = 0; j < k; j++)
            {
                int num = 0;
                double sum = 0;
                for(i = 0; i < len; i++)
                {
                    if(center_pos[i] == j)
                    {
                        num++;
                        sum += can_vec[i];
                    }
                }
                if(num!=0)
                    centers[j] = sum/num;
            }
        }
        for(i = 0; i < k; i++)
            cout << centers[i] << "  ";
        cout << endl;
        for(i = 0; i < len; i++)
        {
            err_dis += fabs(can_vec[i]-centers[center_pos[i]]);
        }
        cout << k << " " << err_dis << endl;
        if(err_dis < min_err_dis)
        {
            min_err_dis = err_dis;
            min_error_k = k;
        }
        for(i = 0; i < k; i++)
        {
            ans[k].push_back(centers[i]);
        }
    }
    cout << "min_error_k: " << min_error_k << endl;
    for(i = 0; i < min_error_k; i++)
    {
        ans_vec.push_back(ans[min_error_k][i]);
    }
    FILE *fp = fopen("mid_num.txt", "w");
    for(i = 0; i < len; i++)
    {
        fprintf(fp, "%d\n", center_pos[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
}
*/

void SLIC::DrawContoursAroundSegments(
    unsigned char*			ubuff,
    const int*				labels,
    const int&				width,
    const int&				height,
    const cv::Scalar&		color )
{
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    int sz = width*height;

    vector<bool> istaken(sz, false);

    int mainindex(0);
    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            int np(0);
            for( int i = 0; i < 8; i++ )
            {
                int x = k + dx8[i];
                int y = j + dy8[i];

                if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                {
                    int index = y*width + x;

                    if( false == istaken[index] )//comment this to obtain internal contours
                    {
                        if( labels[mainindex] != labels[index] ) np++;
                    }
                }
            }
            if( np > 1 )//change to 2 or 3 for thinner lines
            {
                ubuff[mainindex] = (uchar)color.val[0];
                istaken[mainindex] = true;
            }
            mainindex++;
        }
    }
}

//=================================================================================
/// DrawContoursAroundSegmentsTwoColors
///
/// Internal contour drawing option exists. One only needs to comment the if
/// statement inside the loop that looks at neighbourhood.
//=================================================================================
void SLIC::DrawContoursAroundSegmentsTwoColors(
    unsigned int*			img,
    const int*				labels,
    const int&				width,
    const int&				height)
{
    const int dx[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    int sz = width*height;

    vector<bool> istaken(sz, false);

    vector<int> contourx(sz);
    vector<int> contoury(sz);
    int mainindex(0);
    int cind(0);
    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            int np(0);
            for( int i = 0; i < 8; i++ )
            {
                int x = k + dx[i];
                int y = j + dy[i];

                if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                {
                    int index = y*width + x;

                    //if( false == istaken[index] )//comment this to obtain internal contours
                    {
                        if( labels[mainindex] != labels[index] ) np++;
                    }
                }
            }
            if( np > 1 )
            {
                contourx[cind] = k;
                contoury[cind] = j;
                istaken[mainindex] = true;
                //img[mainindex] = color;
                cind++;
            }
            mainindex++;
        }
    }

    int numboundpix = cind;//int(contourx.size());

    for( int j = 0; j < numboundpix; j++ )
    {
        int ii = contoury[j]*width + contourx[j];
        img[ii] = 0xffffff;
        //----------------------------------
        // Uncomment this for thicker lines
        //----------------------------------
        for( int n = 0; n < 8; n++ )
        {
            int x = contourx[j] + dx[n];
            int y = contoury[j] + dy[n];
            if( (x >= 0 && x < width) && (y >= 0 && y < height) )
            {
                int ind = y*width + x;
                if(!istaken[ind]) img[ind] = 0;
            }
        }
    }
}


//==============================================================================
///	DetectLabEdges
//==============================================================================
void SLIC::DetectLabEdges(
    const double*				lvec,
    const double*				avec,
    const double*				bvec,
    const int&					width,
    const int&					height,
    vector<double>&				edges)
{
    int sz = width*height;

    edges.resize(sz,0);
    for( int j = 1; j < height-1; j++ )
    {
        for( int k = 1; k < width-1; k++ )
        {
            int i = j*width+k;

            double dx = (lvec[i-1]-lvec[i+1])*(lvec[i-1]-lvec[i+1]) +
                        (avec[i-1]-avec[i+1])*(avec[i-1]-avec[i+1]) +
                        (bvec[i-1]-bvec[i+1])*(bvec[i-1]-bvec[i+1]);

            double dy = (lvec[i-width]-lvec[i+width])*(lvec[i-width]-lvec[i+width]) +
                        (avec[i-width]-avec[i+width])*(avec[i-width]-avec[i+width]) +
                        (bvec[i-width]-bvec[i+width])*(bvec[i-width]-bvec[i+width]);

            //edges[i] = (sqrt(dx) + sqrt(dy));
            edges[i] = (dx + dy);
        }
    }
}

//===========================================================================
///	PerturbSeeds
//===========================================================================
void SLIC::PerturbSeeds(
    vector<double>&				kseedsl,
    vector<double>&				kseedsa,
    vector<double>&				kseedsb,
    vector<double>&				kseedsx,
    vector<double>&				kseedsy,
    const vector<double>&		edges)
{
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    int numseeds = kseedsl.size();

    for( int n = 0; n < numseeds; n++ )
    {
        int ox = kseedsx[n];//original x
        int oy = kseedsy[n];//original y
        int oind = oy*m_width + ox;

        int storeind = oind;
        for( int i = 0; i < 8; i++ )
        {
            int nx = ox+dx8[i];//new x
            int ny = oy+dy8[i];//new y

            if( nx >= 0 && nx < m_width && ny >= 0 && ny < m_height)
            {
                int nind = ny*m_width + nx;
                if( edges[nind] < edges[storeind])
                {
                    storeind = nind;
                }
            }
        }
        if(storeind != oind)
        {
            kseedsx[n] = storeind%m_width;
            kseedsy[n] = storeind/m_width;
            kseedsl[n] = m_lvec[storeind];
            kseedsa[n] = m_avec[storeind];
            kseedsb[n] = m_bvec[storeind];
        }
    }
}


//===========================================================================
///	GetLABXYSeeds_ForGivenStepSize
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================
void SLIC::GetLABXYSeeds_ForGivenStepSize(
    vector<double>&				kseedsl,
    vector<double>&				kseedsa,
    vector<double>&				kseedsb,
    vector<double>&				kseedsx,
    vector<double>&				kseedsy,
    const int&					STEP,
    const bool&					perturbseeds,
    const vector<double>&		edgemag)
{
    int numseeds(0);
    int n(0);

    //int xstrips = m_width/STEP;
    //int ystrips = m_height/STEP;
    int xstrips = (0.5+double(m_width)/double(STEP));
    int ystrips = (0.5+double(m_height)/double(STEP));

    int xerr = m_width  - STEP*xstrips;
    int yerr = m_height - STEP*ystrips;

    double xerrperstrip = double(xerr)/double(xstrips);
    double yerrperstrip = double(yerr)/double(ystrips);

    int xoff = STEP/2;
    int yoff = STEP/2;
    //-------------------------
    numseeds = xstrips*ystrips;
    //-------------------------
    kseedsl.resize(numseeds);
    kseedsa.resize(numseeds);
    kseedsb.resize(numseeds);
    kseedsx.resize(numseeds);
    kseedsy.resize(numseeds);

    for( int y = 0; y < ystrips; y++ )
    {
        int ye = y*yerrperstrip;
        for( int x = 0; x < xstrips; x++ )
        {
            int xe = x*xerrperstrip;
            int i = (y*STEP+yoff+ye)*m_width + (x*STEP+xoff+xe);

            kseedsl[n] = m_lvec[i];
            kseedsa[n] = m_avec[i];
            kseedsb[n] = m_bvec[i];
            kseedsx[n] = (x*STEP+xoff+xe);
            kseedsy[n] = (y*STEP+yoff+ye);
            n++;
        }
    }


    if(perturbseeds)
    {
        PerturbSeeds(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, edgemag);
    }
}

//===========================================================================
///	GetLABXYSeeds_ForGivenK
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================
void SLIC::GetLABXYSeeds_ForGivenK(
    vector<double>&				kseedsl,
    vector<double>&				kseedsa,
    vector<double>&				kseedsb,
    vector<double>&				kseedsx,
    vector<double>&				kseedsy,
    const int&					K,
    const bool&					perturbseeds,
    const vector<double>&		edgemag)
{
    int sz = m_width*m_height;
    double step = sqrt(double(sz)/double(K));
    int T = step;
    int xoff = step/2;
    int yoff = step/2;

    int n(0);
    int r(0);
    for( int y = 0; y < m_height; y++ )
    {
        int Y = y*step + yoff;
        if( Y > m_height-1 ) break;

        for( int x = 0; x < m_width; x++ )
        {
            //int X = x*step + xoff;//square grid
            int X = x*step + (xoff<<(r&0x1));//hex grid
            if(X > m_width-1) break;

            int i = Y*m_width + X;

            //_ASSERT(n < K);

            //kseedsl[n] = m_lvec[i];
            //kseedsa[n] = m_avec[i];
            //kseedsb[n] = m_bvec[i];
            //kseedsx[n] = X;
            //kseedsy[n] = Y;
            kseedsl.push_back(m_lvec[i]);
            kseedsa.push_back(m_avec[i]);
            kseedsb.push_back(m_bvec[i]);
            kseedsx.push_back(X);
            kseedsy.push_back(Y);
            n++;
        }
        r++;
    }

    if(perturbseeds)
    {
        PerturbSeeds(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, edgemag);
    }
}


//===========================================================================
///	PerformSuperpixelSegmentation_VariableSandM
///
///	Magic SLIC - no parameters
///
///	Performs k mean segmentation. It is fast because it looks locally, not
/// over the entire image.
/// This function picks the maximum value of color distance as compact factor
/// M and maximum pixel distance as grid step size S from each cluster (13 April 2011).
/// So no need to input a constant value of M and S. There are two clear
/// advantages:
///
/// [1] The algorithm now better handles both textured and non-textured regions
/// [2] There is not need to set any parameters!!!
///
/// SLICO (or SLIC Zero) dynamically varies only the compactness factor S,
/// not the step size S.
//===========================================================================
void SLIC::PerformSuperpixelSegmentation_VariableSandM(
    vector<double>&				kseedsl,
    vector<double>&				kseedsa,
    vector<double>&				kseedsb,
    vector<double>&				kseedsx,
    vector<double>&				kseedsy,
    int*						klabels,
    const int&					STEP,
    const int&					NUMITR)
{
    int sz = m_width*m_height;
    const int numk = kseedsl.size();
    //double cumerr(99999.9);
    int numitr(0);

    //----------------
    int offset = STEP;
    if(STEP < 10) offset = STEP*1.5;
    //----------------

    vector<double> sigmal(numk, 0);
    vector<double> sigmaa(numk, 0);
    vector<double> sigmab(numk, 0);
    vector<double> sigmax(numk, 0);
    vector<double> sigmay(numk, 0);
    vector<int> clustersize(numk, 0);
    vector<double> inv(numk, 0);//to store 1/clustersize[k] values
    vector<double> distxy(sz, DBL_MAX);
    vector<double> distlab(sz, DBL_MAX);
    vector<double> distvec(sz, DBL_MAX);
    vector<double> maxlab(numk, 10*10);//THIS IS THE VARIABLE VALUE OF M, just start with 10
    vector<double> maxxy(numk, STEP*STEP);//THIS IS THE VARIABLE VALUE OF M, just start with 10

    double invxywt = 1.0/(STEP*STEP);//NOTE: this is different from how usual SLIC/LKM works

    while( numitr < NUMITR )
    {
        //------
        //cumerr = 0;
        numitr++;
        //------

        distvec.assign(sz, DBL_MAX);
        for( int n = 0; n < numk; n++ )
        {
            int y1 = std::max(0,		(int)(kseedsy[n]-offset));
            int y2 = std::min(m_height,	(int)(kseedsy[n]+offset));
            int x1 = std::max(0,		(int)(kseedsx[n]-offset));
            int x2 = std::min(m_width,	(int)(kseedsx[n]+offset));

            for( int y = y1; y < y2; y++ )
            {
                for( int x = x1; x < x2; x++ )
                {
                    int i = y*m_width + x;
                    //_ASSERT( y < m_height && x < m_width && y >= 0 && x >= 0 );
                    assert(y < m_height && x < m_width && y >= 0 && x >= 0);

                    double l = m_lvec[i];
                    double a = m_avec[i];
                    double b = m_bvec[i];

                    distlab[i] =	(l - kseedsl[n])*(l - kseedsl[n]) +
                                    (a - kseedsa[n])*(a - kseedsa[n]) +
                                    (b - kseedsb[n])*(b - kseedsb[n]);

                    distxy[i] =		(x - kseedsx[n])*(x - kseedsx[n]) +
                                    (y - kseedsy[n])*(y - kseedsy[n]);

                    //------------------------------------------------------------------------
                    double dist = distlab[i]/maxlab[n] + distxy[i]*invxywt;//only varying m, prettier superpixels
                    //double dist = distlab[i]/maxlab[n] + distxy[i]/maxxy[n];//varying both m and S
                    //------------------------------------------------------------------------

                    if( dist < distvec[i] )
                    {
                        distvec[i] = dist;
                        klabels[i]  = n;
                    }
                }
            }
        }
        //-----------------------------------------------------------------
        // Assign the max color distance for a cluster
        //-----------------------------------------------------------------
        if(0 == numitr)
        {
            maxlab.assign(numk,1);
            maxxy.assign(numk,1);
        }
        {
            for( int i = 0; i < sz; i++ )
            {
                if(maxlab[klabels[i]] < distlab[i]) maxlab[klabels[i]] = distlab[i];
                if(maxxy[klabels[i]] < distxy[i]) maxxy[klabels[i]] = distxy[i];
            }
        }
        //-----------------------------------------------------------------
        // Recalculate the centroid and store in the seed values
        //-----------------------------------------------------------------
        sigmal.assign(numk, 0);
        sigmaa.assign(numk, 0);
        sigmab.assign(numk, 0);
        sigmax.assign(numk, 0);
        sigmay.assign(numk, 0);
        clustersize.assign(numk, 0);

        for( int j = 0; j < sz; j++ )
        {
            int temp = klabels[j];
            //_ASSERT(klabels[j] >= 0);
            assert(klabels[j] >= 0);
            sigmal[klabels[j]] += m_lvec[j];
            sigmaa[klabels[j]] += m_avec[j];
            sigmab[klabels[j]] += m_bvec[j];
            sigmax[klabels[j]] += (j%m_width);
            sigmay[klabels[j]] += (j/m_width);

            clustersize[klabels[j]]++;
        }

        {
            for( int k = 0; k < numk; k++ )
            {
                //_ASSERT(clustersize[k] > 0);
                if( clustersize[k] <= 0 ) clustersize[k] = 1;
                inv[k] = 1.0/double(clustersize[k]);//computing inverse now to multiply, than divide later
            }
        }

        {
            for( int k = 0; k < numk; k++ )
            {
                kseedsl[k] = sigmal[k]*inv[k];
                kseedsa[k] = sigmaa[k]*inv[k];
                kseedsb[k] = sigmab[k]*inv[k];
                kseedsx[k] = sigmax[k]*inv[k];
                kseedsy[k] = sigmay[k]*inv[k];
            }
        }
    }
}


void SLIC::GetRegionInfo(int x, int y)
{
    int i, j;
    int index = y*m_width+x;
    int my_label = label[index];
    //cout << "my_label: " << my_label << endl;
    int sum = 0;
    double suml = 0, suma = 0, sumb = 0;
    for(i = 0; i < m_height; i++)
    {
        for(j = 0; j < m_width; j++)
        {
            index = i*m_width+j;
            if(label[index] == my_label)
            {
                suml += m_lvec[index];
                suma += m_avec[index];
                sumb += m_bvec[index];
                sum++;
            }
        }
    }
    cout << suml/sum << " " << suma/sum << " " << sumb/sum << endl;
}

//===========================================================================
///	SaveSuperpixelLabels
///
///	Save labels in raster scan order.
//===========================================================================
void SLIC::SaveSuperpixelLabels(
    const int*					labels,
    const int&					width,
    const int&					height,
    const string&				filename,
    const string&				path)
{
    int sz = width*height;
    /*
    int _MAX_FNAME = 1000;

    char fname[_MAX_FNAME];
    char extn[_MAX_FNAME];
    _splitpath(filename.c_str(), NULL, NULL, fname, extn);
    */
    string temp = "label_data";//fname;

    ofstream outfile;
    string finalpath = path + temp + string(".dat");
    outfile.open(finalpath.c_str(), ios::binary);
    for( int i = 0; i < sz; i++ )
    {
        outfile.write((const char*)&labels[i], sizeof(int));
    }
    outfile.close();
}

//===========================================================================
///	EnforceLabelConnectivity
///
///		1. finding an adjacent label for each new component at the start
///		2. if a certain component is too small, assigning the previously found
///		    adjacent label to this component, and not incrementing the label.
//===========================================================================
void SLIC::EnforceLabelConnectivity(
    const int*					labels,//input labels that need to be corrected to remove stray labels
    const int&					width,
    const int&					height,
    int*						nlabels,//new labels
    int&						numlabels,//the number of labels changes in the end if segments are removed
    const int&					K) //the number of superpixels desired by the user
{
    //	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    //	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    const int dx4[4] = {-1,  0,  1,  0};
    const int dy4[4] = { 0, -1,  0,  1};

    const int sz = width*height;
    const int SUPSZ = sz/K;
    //nlabels.resize(sz, -1);
    for( int i = 0; i < sz; i++ ) nlabels[i] = -1;
    int label(0);
    int* xvec = new int[sz];
    int* yvec = new int[sz];
    int oindex(0);
    int adjlabel(0);//adjacent label
    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            if( 0 > nlabels[oindex] )
            {
                nlabels[oindex] = label;
                //--------------------
                // Start a new segment
                //--------------------
                xvec[0] = k;
                yvec[0] = j;
                //-------------------------------------------------------
                // Quickly find an adjacent label for use later if needed
                //-------------------------------------------------------
                {
                    for( int n = 0; n < 4; n++ )
                    {
                        int x = xvec[0] + dx4[n];
                        int y = yvec[0] + dy4[n];
                        if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                        {
                            int nindex = y*width + x;
                            if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
                        }
                    }
                }

                int count(1);
                for( int c = 0; c < count; c++ )
                {
                    for( int n = 0; n < 4; n++ )
                    {
                        int x = xvec[c] + dx4[n];
                        int y = yvec[c] + dy4[n];

                        if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                        {
                            int nindex = y*width + x;

                            if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
                            {
                                xvec[count] = x;
                                yvec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }

                    }
                }
                //-------------------------------------------------------
                // If segment size is less then a limit, assign an
                // adjacent label found before, and decrement label count.
                //-------------------------------------------------------
                if(count <= SUPSZ >> 2)
                {
                    for( int c = 0; c < count; c++ )
                    {
                        int ind = yvec[c]*width+xvec[c];
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }
    numlabels = label;

    if(xvec) delete [] xvec;
    if(yvec) delete [] yvec;
}

//===========================================================================
///	PerformSLICO_ForGivenStepSize
///
/// There is option to save the labels if needed.
//===========================================================================
void SLIC::PerformSLICO_ForGivenStepSize(
    const unsigned int*			ubuff,
    const int					width,
    const int					height,
    int*						klabels,
    int&						numlabels,
    const int&					STEP,
    const double&				m)
{
    vector<double> kseedsl(0);
    vector<double> kseedsa(0);
    vector<double> kseedsb(0);
    vector<double> kseedsx(0);
    vector<double> kseedsy(0);

    //--------------------------------------------------
    m_width  = width;
    m_height = height;
    int sz = m_width*m_height;
    //klabels.resize( sz, -1 );
    //--------------------------------------------------
    //klabels = new int[sz];
    for( int s = 0; s < sz; s++ ) klabels[s] = -1;
    //--------------------------------------------------
    DoRGBtoLABConversion(ubuff, m_lvec, m_avec, m_bvec);
    //--------------------------------------------------

    bool perturbseeds(true);
    vector<double> edgemag(0);
    if(perturbseeds) DetectLabEdges(m_lvec, m_avec, m_bvec, m_width, m_height, edgemag);
    GetLABXYSeeds_ForGivenStepSize(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, STEP, perturbseeds, edgemag);

    PerformSuperpixelSegmentation_VariableSandM(kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,klabels,STEP,10);
    numlabels = kseedsl.size();

    int* nlabels = new int[sz];
    EnforceLabelConnectivity(klabels, m_width, m_height, nlabels, numlabels, double(sz)/double(STEP*STEP));
    {
        for(int i = 0; i < sz; i++ ) klabels[i] = nlabels[i];
    }
    if(nlabels) delete [] nlabels;
}


//===========================================================================
///	PerformSLICO_ForGivenK
///
/// Zero parameter SLIC algorithm for a given number K of superpixels.
//===========================================================================
void SLIC::PerformSLICO_ForGivenK(
    const unsigned int*			ubuff,
    const int					width,
    const int					height,
    int*						klabels,
    int&						numlabels,
    const int&					K,//required number of superpixels
    const double&				m)//weight given to spatial distance
{
    vector<double> kseedsl(0);
    vector<double> kseedsa(0);
    vector<double> kseedsb(0);
    vector<double> kseedsx(0);
    vector<double> kseedsy(0);

    //--------------------------------------------------
    m_width  = width;
    m_height = height;
    int sz = m_width*m_height;
    //--------------------------------------------------
    //if(0 == klabels) klabels = new int[sz];
    for( int s = 0; s < sz; s++ ) klabels[s] = -1;
    //--------------------------------------------------
    if(1)//LAB
    {
        DoRGBtoLABConversion(ubuff, m_lvec, m_avec, m_bvec);
    }
    else//RGB
    {
        m_lvec = new double[sz];
        m_avec = new double[sz];
        m_bvec = new double[sz];
        for( int i = 0; i < sz; i++ )
        {
            m_lvec[i] = ubuff[i] >> 16 & 0xff;
            m_avec[i] = ubuff[i] >>  8 & 0xff;
            m_bvec[i] = ubuff[i]       & 0xff;
        }
    }
    //--------------------------------------------------

    bool perturbseeds(true);
    vector<double> edgemag(0);
    if(perturbseeds) DetectLabEdges(m_lvec, m_avec, m_bvec, m_width, m_height, edgemag);
    GetLABXYSeeds_ForGivenK(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, K, perturbseeds, edgemag);

    int STEP = sqrt(double(sz)/double(K)) + 2.0;//adding a small value in the even the STEP size is too small.
    //PerformSuperpixelSLIC(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, klabels, STEP, edgemag, m);
    PerformSuperpixelSegmentation_VariableSandM(kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,klabels,STEP,10);
    numlabels = kseedsl.size();

    int* nlabels = new int[sz];
    EnforceLabelConnectivity(klabels, m_width, m_height, nlabels, numlabels, K);
    {
        for(int i = 0; i < sz; i++ ) klabels[i] = nlabels[i];
    }
    if(nlabels) delete [] nlabels;
}

void SLIC::PerformSLICO_ForGivenK(
    const unsigned char*		ubuff,
    const int					width,
    const int					height,
    int*						klabels,
    int&						numlabels,
    const int&					K,//required number of superpixels
    const double&				m)//weight given to spatial distance
{
    vector<double> kseedsl(0);
    vector<double> kseedsa(0);
    vector<double> kseedsb(0);
    vector<double> kseedsx(0);
    vector<double> kseedsy(0);

    //--------------------------------------------------
    m_width  = width;
    m_height = height;
    int sz = m_width*m_height;
    //--------------------------------------------------
    //if(0 == klabels) klabels = new int[sz];
    for( int s = 0; s < sz; s++ ) klabels[s] = -1;
    //--------------------------------------------------

    m_lvec = new double[sz];
    m_avec = new double[sz];
    m_bvec = new double[sz];
    for( int i = 0; i < sz; i++ )
    {
        m_lvec[i] = ubuff[i];
        m_avec[i] = 0;
        m_bvec[i] = 0;
    }


    //--------------------------------------------------

    bool perturbseeds(true);
    vector<double> edgemag(0);
    if(perturbseeds) DetectLabEdges(m_lvec, m_avec, m_bvec, m_width, m_height, edgemag);
    GetLABXYSeeds_ForGivenK(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, K, perturbseeds, edgemag);

    int STEP = sqrt(double(sz)/double(K)) + 2.0;//adding a small value in the even the STEP size is too small.
    //PerformSuperpixelSLIC(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, klabels, STEP, edgemag, m);
    PerformSuperpixelSegmentation_VariableSandM(kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,klabels,STEP,10);
    numlabels = kseedsl.size();

    int* nlabels = new int[sz];
    EnforceLabelConnectivity(klabels, m_width, m_height, nlabels, numlabels, K);
    {
        for(int i = 0; i < sz; i++ ) klabels[i] = nlabels[i];
    }
    if(nlabels) delete [] nlabels;
}

void SLIC::GenerateSuperpixels(cv::Mat& img, UINT numSuperpixels)
{
    if (img.empty())
    {
        exit(-1);
    }

    int height = img.rows;
    int width = img.cols;
    int sz = height * width;

    label = new int[sz];
    point_flag = new int[sz];
    edge_flag = new int[sz];
    m_supers = numSuperpixels;
    //cout << "channels: " << img.channels() << endl;
    if (img.channels() == 1)
    {
        type = GRAY;
    }
    else if (img.channels() == 3)
    {
        type = RGB;
    }
    if (type == GRAY)
    {
        Mat2Buffer(img, bufferGray);
        PerformSLICO_ForGivenK(bufferGray, img.cols, img.rows, label, sz, numSuperpixels, 10);
    }
    else if (type == RGB)
    {
        Mat2Buffer(img, bufferRGB);
        PerformSLICO_ForGivenK(bufferRGB, img.cols, img.rows, label, sz, numSuperpixels, 10);
    }
}

//
int* SLIC::GetLabel()
{
    return label;
}

void SLIC::output_label()
{
    FILE *fp = fopen("labels.txt", "w");
    int x, y;
    //cout << m_height << " " << m_width << endl;
    for(y = 0; y < m_height; y++)
    {
        for(x = 0; x < m_width; x++)
        {
            fprintf(fp, "%d ", label[y*m_width+x]);
        }
        fprintf(fp, "\n");
    }
}
void SLIC::GenerateAveValue()
{
    int i,j,k;

    for(k = 0; k < m_supers; k++)
    {
        int sum = 0;
        double suml = 0, suma = 0, sumb = 0;
        for(i = 0; i < m_height; i++)
        {
            for(j = 0; j < m_width; j++)
            {
                int index = i*m_width+j;
                if(label[index] == k)
                {
                    suml += m_lvec[index];
                    suma += m_avec[index];
                    sumb += m_bvec[index];
                    sum++;
                }
            }
        }
        m_lave[k] = suml/sum;
        m_aave[k] = suma/sum;
        m_bave[k] = sumb/sum;
    }
}
//根据分割边界值，对每一块进行标记
void SLIC::MarkLabel()
{
    int k;
    //m_partition = 12;
    for(k = 0; k < m_supers; k++)
    {
        if(m_aave[k] < m_partition)
            label_flag[k] = 1;
        else
            label_flag[k] = 0;
    }
}
//寻找分割的边界值
//对lower和upper的区域内超像素块的值做一个聚类，分为2类，第一类为接近牙齿区域的a值，第二类为接近牙龈区域的a值
int SLIC::GeneratePartitionValue()
{
    int i, j, k;
    int candidates_flag[1000];
    memset(candidates_flag, 0, sizeof(candidates_flag));
    int lower, upper, left_, right_;
    GetLowerUpper(m_height, m_width, lower, upper);
    //lower = m_height/5, upper = 4*lower;
    left_ = 250, right_ = 0.8*m_width;
    for(i = lower; i <= upper; i++)
    {
        for(j = left_; j <= right_; j++)
        {
            int index = i*m_width+j;
            int my_label = label[index];
            candidates_flag[my_label] = 1;
        }
    }
    int par_pos[1000];
    memset(par_pos, 0, sizeof(par_pos));
    j = 0;
    for(i = 0; i < m_supers; i++)
    {
        if(candidates_flag[i] == 1)
            m_par[j++] = m_aave[i];
    }
    //cout << j << endl;
    sort(m_par, m_par+j);

    double center1 = m_par[0], center2 = m_par[j-1];
    int ite = 0;
    //cout << "clustering...\n";
    while(!ite)
    {
        double dist1, dist2;
        int new_pos;
        ite = 1;
        for(i = 0; i < j; i++)
        {
            dist1 = fabs(m_par[i]-center1);
            dist2 = fabs(m_par[i]-center2);
            new_pos = dist1<dist2? 1:2;
            if(new_pos!=par_pos[i])
                ite = 0;
            par_pos[i] = new_pos;
        }
        double sum1 = 0.0, sum2 = 0.0;
        int num1 = 0.0, num2 = 0;
        for(i = 0; i < j; i++)
        {
            if(par_pos[i] == 1)
            {
                num1++;
                sum1+=m_par[i];
            }
            else
            {
                num2++;
                sum2+=m_par[i];
            }
        }
        center1 = sum1/num1;
        center2 = sum2/num2;
    }
    FILE *fp = fopen("m_aave.txt", "w");
    for(i = 0; i < j; i++)
    {
        fprintf(fp, "%lf %d\n", m_par[i], par_pos[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    m_partition = m_par[j-1];
    for(i = 0; i < j-1; i++)
    {
        if(par_pos[i]==1 && par_pos[i+1] == 2)
        {
            m_partition = m_par[i];
            break;
        }
    }
    //m_partition = 18;
    //Debug(m_partition);
    if(i==0)
        i = j-1;
    return i;

}
int SLIC::GetImgWithContours(cv::Scalar color)
{
    //cout << type << " " << GRAY << endl;
    //output_label();
    int i, j, k;
    int ans[4];
    int prior_ans = 0;
    memset(ans, 0, sizeof(ans));
    GenerateAveValue();
    int ite_sum = GeneratePartitionValue();
    //Debug(ite_sum);
    int times = 0;
    for(i = ite_sum-1; i >= 2*ite_sum/3; i--)
    {//cout << i << endl;
        times++;
        if(times >= 50)
            break;
        m_partition = m_par[i];
        MarkLabel();
        if (type == GRAY)
        {
            DrawContoursAroundSegments(bufferGray, label, m_width, m_height, color);
            /*
            cv::Mat result(m_height, m_width, CV_8UC1);
            memcpy(result.data, bufferGray, m_width*m_height*sizeof(uchar));
            return result;
            */
        }
        else if (type == RGB)
        {
            DrawContoursAroundSegments(bufferRGB, label, m_width, m_height, color);
            /*
            cv::Mat result(m_height, m_width, CV_8UC4);
            memcpy(result.data, bufferRGB, m_width*m_height*sizeof(UINT));
            cvtColor(result, result, CV_BGRA2BGR);
            return result;
            */
        }
        int ans_num = TestEdge(m_height, m_width);
        ans[ans_num]++;
        if(i == ite_sum-1)
        {
            prior_ans = ans_num;
        }
    }
    int max_num;// = max(ans[1], max(ans[2], ans[3]));
    
    if(ans[1] > ans[2])
    {
        if(ans[1] > ans[3])
            max_num = 1;
        else if(ans[1] < ans[3])
            max_num = 3;
        else
            max_num = prior_ans;
    }
    else if(ans[1] < ans[2])
    {
        if(ans[2] > ans[3])
            max_num = 2;
        else if(ans[2] < ans[3])
            max_num = 3;
        else
            max_num = prior_ans;
    }
    else if(ans[1] < ans[3])
        max_num = 3;
    else
        max_num = prior_ans;
    cout << max_num << endl;
    return max_num;

}

void SLIC::Mat2Buffer(const cv::Mat& img, UINT*& buffer)
{
    int sz = img.cols * img.rows;
    if (bufferRGB)
    {
        delete [] bufferRGB;
    }
    bufferRGB    = new UINT[sz];

    cv::Mat newImage;
    cv::cvtColor(img, newImage, CV_BGR2BGRA);
    memcpy( bufferRGB, (UINT*)newImage.data, sz*sizeof(UINT) );

}

void SLIC::Mat2Buffer(const cv::Mat& img, uchar*& buffer)
{
    int sz = img.cols * img.rows;
    if (bufferGray)
    {
        delete [] bufferGray;
    }
    bufferGray    = new uchar[sz];

    memcpy( bufferGray, (UINT*)img.data, sz*sizeof(uchar) );

}
