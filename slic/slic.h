// SLIC.h: interface for the SLIC class.
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

#if !defined(_SLIC_H_INCLUDED_)
#define _SLIC_H_INCLUDED_


#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace cv;

#define Debug(x) cout<<#x<<"="<<(x)<<endl;

typedef unsigned int UINT;
typedef unsigned char uchar;

enum imageType {RGB, GRAY};

class SLIC  
{
public:
	SLIC();
	virtual ~SLIC();

	//===========================================================================
	///	Perform SLIC algorithm on the given image 
	/// with the given number of superpixels
	//===========================================================================
	void GenerateSuperpixels(
		cv::Mat& img, 
		UINT numSuperpixels);

	//===========================================================================
	///	Get label on each pixel which shows the number of superpixel it belongs to
	//===========================================================================
	int* GetLabel();

	//===========================================================================
	///	Get the result image with contours on the given color
	//===========================================================================
	int GetImgWithContours(cv::Scalar color);

	void GetRegionInfo(int x, int y);
	void MarkLabel();
	int TestEdge(int height, int width);

private:

	//============================================================================
	// Superpixel segmentation for a given step size (superpixel size ~= step*step)
	//============================================================================
	void PerformSLICO_ForGivenStepSize(
		const unsigned int*			ubuff,//Each 32 bit unsigned int contains ARGB pixel values.
		const int					width,
		const int					height,
		int*						klabels,
		int&						numlabels,
		const int&					STEP,
		const double&				m);
	//============================================================================
	// Superpixel segmentation for a given number of superpixels
	//============================================================================
	void PerformSLICO_ForGivenK(
		const unsigned int*			ubuff,//Each 32 bit unsigned int contains ARGB pixel values.
		const int					width,
		const int					height,
		int*						klabels,
		int&						numlabels,
		const int&					K,
		const double&				m);

	void PerformSLICO_ForGivenK(
		const unsigned char*		ubuff,
		const int					width,
		const int					height,
		int*						klabels,
		int&						numlabels,
		const int&					K,//required number of superpixels
		const double&				m);//weight given to spatial distance

	//============================================================================
	// Save superpixel labels in a text file in raster scan order
	//============================================================================
	void SaveSuperpixelLabels(
		const int*					labels,
		const int&					width,
		const int&					height,
		const string&				filename,
		const string&				path);
	//============================================================================
	// Function to draw boundaries around superpixels of a given 'color'.
	// Can also be used to draw boundaries around supervoxels, i.e layer by layer.
	//============================================================================
	void DrawContoursAroundSegments(
		unsigned int*				ubuff,
		const int*					labels,
		const int&					width,
		const int&					height,
		const cv::Scalar&			color );

	void DrawContoursAroundSegments(
		unsigned char*			ubuff,
		const int*				labels,
		const int&				width,
		const int&				height,
		const cv::Scalar&		color );

	void DrawContoursAroundSegmentsTwoColors(
		unsigned int*				ubuff,
		const int*					labels,
		const int&					width,
		const int&					height);




private:

	//============================================================================
	// Magic SLIC. No need to set M (compactness factor) and S (step size).
	// SLICO (SLIC Zero) varies only M dynamicaly, not S.
	//============================================================================
	void PerformSuperpixelSegmentation_VariableSandM(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		int*						klabels,
		const int&					STEP,
		const int&					NUMITR);
	//============================================================================
	// Pick seeds for superpixels when step size of superpixels is given.
	//============================================================================
	void GetLABXYSeeds_ForGivenStepSize(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&		edgemag);
	//============================================================================
	// Pick seeds for superpixels when number of superpixels is input.
	//============================================================================
	void GetLABXYSeeds_ForGivenK(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&		edges);

	//============================================================================
	// Move the seeds to low gradient positions to avoid putting seeds at region boundaries.
	//============================================================================
	void PerturbSeeds(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const vector<double>&		edges);
	//============================================================================
	// Detect color edges, to help PerturbSeeds()
	//============================================================================
	void DetectLabEdges(
		const double*				lvec,
		const double*				avec,
		const double*				bvec,
		const int&					width,
		const int&					height,
		vector<double>&				edges);
	//============================================================================
	// xRGB to XYZ conversion; helper for RGB2LAB()
	//============================================================================
	void RGB2XYZ(
		const int&					sR,
		const int&					sG,
		const int&					sB,
		double&						X,
		double&						Y,
		double&						Z);
	//============================================================================
	// sRGB to CIELAB conversion
	//============================================================================
	void RGB2LAB(
		const int&					sR,
		const int&					sG,
		const int&					sB,
		double&						lval,
		double&						aval,
		double&						bval);
	//============================================================================
	// sRGB to CIELAB conversion for 2-D images
	//============================================================================
	void DoRGBtoLABConversion(
		const unsigned int*&		ubuff,
		double*&					lvec,
		double*&					avec,
		double*&					bvec);
	//============================================================================
	// sRGB to CIELAB conversion for 3-D volumes
	//============================================================================
	void DoRGBtoLABConversion(
		const unsigned int**&		ubuff,
		double**&					lvec,
		double**&					avec,
		double**&					bvec);

	//============================================================================
	// Post-processing of SLIC segmentation, to avoid stray labels.
	//============================================================================
	void EnforceLabelConnectivity(
		const int*					labels,
		const int&					width,
		const int&					height,
		int*						nlabels,//input labels that need to be corrected to remove stray labels
		int&						numlabels,//the number of labels changes in the end if segments are removed
		const int&					K); //the number of superpixels desired by the user

	void Mat2Buffer(const cv::Mat& img, UINT*& buffer);

	void Mat2Buffer(const cv::Mat& img, uchar*& buffer);
	void output_label();
	int GeneratePartitionValue();
	void GenerateAveValue();
	void SmartCluster(vector<double> can_vec, vector<double> &ans_vec);
	void GetLowerUpper(int height, int width, int &lower, int &upper);
	int GetClassification(vector<int> ans_up, vector<int> ans_down);




private:
	int										m_width;
	int										m_height;
	int										m_depth;

	int 									m_supers;

	double*									m_lvec;
	double*									m_avec;
	double*									m_bvec;

	double**								m_lvecvec;
	double**								m_avecvec;
	double**								m_bvecvec;

	UINT*									bufferRGB; // buffer for if RGB image
	uchar*									bufferGray; // buffer if gray image

	int*									label; // label record the superpixel pixel belongs to
	int*									point_flag;//0 red; 1 green
	int* 									edge_flag;//0 not edge; 1 edge
	imageType								type;

	int 									label_flag[5000];
	double 									m_lave[5000], m_aave[5000], m_bave[5000];
	double									m_partition;
	double 									m_par[1000];
	vector<int> 							mid_vec_up[200], mid_vec_down[200];
	//vector<int> 							mid_y_up[10], mid_y_down[10];


};

#endif // !defined(_SLIC_H_INCLUDED_)
