#include <iostream>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <list>
#include <vector>
#include <cstdlib>
#include <forward_list>
#include <cmath>
#include <time.h>
#include "highgui.h"
#include "cv.h"

#include "Slic.h"
#include "SlicMerge.h"
#include "HierarchicalTree.h"

using namespace cv;

#pragma comment(linker, "/STACK:102400000,102400000")    //防止栈溢出

int width;
int height;

int main()
{

    int sz;
    int i, ii;
    int x, y;

    int step;

    int numseeds;

    int k;

    int dims[2] = {0};
    int* outputNumSuperpixels;
    int* outlabels;
    int finalNumberOfLabels;
    unsigned char* imgbytes; //等同于src.data
    int numelements;
    int numSuperpixels = 200;//default value
    double compactness = 10;//default value
	double maxDiffence = 200;//default value
	
	numSuperpixels = 3000; //**超像素个数,适用于demo
    compactness = 10; //**紧凑度
	maxDiffence = 20; //**允许的最大异质性数值

	Mat zy1, zy2, zy3, zy4, srimg;
	zy1 = imread("C:/b.bmp",0);
	zy2 = imread("C:/g.bmp",0);
	zy3 = imread("C:/r.bmp",0);
	zy4 = imread("C:/ir.bmp",0);
	if (zy1.empty() || zy2.empty() || zy3.empty() || zy4.empty())
	{
		printf("Can not open Image\n");
		system("pause");
		exit(0);
	}

	CvSize size;
	size.width = zy1.cols;
	size.height = zy1.rows;
	srimg.create(size.height, size.width, CV_8UC4);
	for (int i = 0;i<size.height;i++)
		for (int j = 0; j<size.width; j++)
		{
			srimg.data[(i*size.width+j)*4]=zy3.data[i*size.width+j];
			srimg.data[(i*size.width+j)*4+1]=zy2.data[i*size.width+j];
			srimg.data[(i*size.width+j)*4+2]=zy1.data[i*size.width+j];
			srimg.data[(i*size.width+j)*4+3]=zy4.data[i*size.width+j];
		}
	namedWindow("srimg");
	imshow("srimg",srimg);
	waitKey(0);
    //分离三个RGB通道
    width = srimg.cols;
    height = srimg.rows;
    sz = width*height;
    int* rin = new int[height*width];       //@可以用uchar节省空间
    int* gin = new int[height*width];
    int* bin = new int[height*width];
	int* nirIn = new int[height*width];		//近红外波段
    double* lvec = new double[height*width]; //@同时为LAB开辟空间，可以优化
    double* avec = new double[height*width];
    double* bvec = new double[height*width];
    int* klabels = new int[height*width];
    int* clabels = new int[height*width];
    int* seedIndices = new int[height*width];


	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
		{
			bin[i*width + j] = srimg.data[(i*width + j)*4];
			gin[i*width + j] = srimg.data[(i*width + j)*4 + 1];
			rin[i*width + j] = srimg.data[(i*width + j)*4 + 2];
			nirIn[i*width+j] = srimg.data[(i*width + j)*4 + 3];
		}
    //部分参数定义  转换可能有所不同
    numelements = 3;   //默认只对彩色图像做处理
    //numdims明显的二维
    dims[0] = height;
    dims[1] = width;
    imgbytes = srimg.data; //指向原图像数据域

    

    //convert from rgb to lab
    rgbtolab(rin,gin,bin,sz,lvec,avec,bvec);
    
    //find seeds
    step = sqrt((double)(sz)/(double)(numSuperpixels))+0.5;
    getLABXYSeeds(step,width,height,seedIndices,&numseeds);

    double* kseedsx = new double[numseeds];
    double* kseedsy = new double[numseeds];
    double* kseedsl = new double[numseeds];
    double* kseedsa = new double[numseeds];
    double* kseedsb = new double[numseeds];
    for(k = 0; k < numseeds; k++)
    {
        kseedsx[k] = seedIndices[k]%width;
        kseedsy[k] = seedIndices[k]/width;
        kseedsl[k] = lvec[seedIndices[k]];
        kseedsa[k] = avec[seedIndices[k]];
        kseedsb[k] = bvec[seedIndices[k]];
    }

    //Compute superpixels

    PerformSuperpixelSLIC(lvec, avec, bvec, kseedsl,kseedsa,kseedsb,kseedsx,kseedsy,width,height,numseeds,klabels,step,compactness);

    //Enforce connectivity
    EnforceSuperpixelConnectivity(klabels,width,height,numSuperpixels,clabels,&finalNumberOfLabels);
    
    //output
    //clabels为连通矩阵 finalNumberOfLabels为超像素区域个数
	printf("%d\n", finalNumberOfLabels);
	for (int i = 0; i<20; i++)
	{
		for (int j = 0; j<9; j++)
		{
			printf("%d ", clabels[i*width + j]);
		}
		printf("\n");
	}
	
	//实例化基层超像素对象
	CSuperPixelSet* csps = new CSuperPixelSet[finalNumberOfLabels]; 
	//实例化层次树各个节点
	BTreeNode* hierarchicalTree = new BTreeNode[2*finalNumberOfLabels-1]; 
	//超像素过分割，建立超像素实体集
	createSuperPixelVector(clabels, width, height, csps, srimg, finalNumberOfLabels, hierarchicalTree);
	//实例化拓扑图头结点
	ArrayHeadGraphNode *mAhgn = new ArrayHeadGraphNode[2*finalNumberOfLabels-1];   //严格的2n-1，最终层次树的头结点下标为2n-2
	//建立初始拓扑图
	createToplogicalGraph(clabels, width, height, mAhgn,numSuperpixels);


	printf("\n最终层次树结点数：%d\n", 2*finalNumberOfLabels - 2);
	system("pause");
	
	/****************************************/
	//构建层次树
	createHierarchicalTree(mAhgn, hierarchicalTree, srimg, maxDiffence, finalNumberOfLabels);
	/***************************************/

	//后处理
	int level;
	int *newLabels = new int[height*width] ();
	do{
		printf("\n\n请输入查询结点等级:");
		scanf("%d", &level);
		//searchTreeNodeWithLevel(&hierarchicalTree[2*finalNumberOfLabels-2], level, finalNumberOfLabels);
		int setValue = -1;
		setAllNodeValue(newLabels, level, &hierarchicalTree[2*finalNumberOfLabels-2], setValue, csps);
		


		printf("objectNum: %d\n", setValue+1);
		
		//放弃层次树结点中的其它信息，只保留层次信息
		//基于新联通图层建立新对象面块集合以及拓扑关系信息
		int objectNum = setValue + 1;
		ObjectNode* oNode = new ObjectNode[objectNum];
		ArrayHeadGraphNode *newAHGn = new ArrayHeadGraphNode[objectNum];
		createNewObjectSet(newLabels, srimg, oNode, objectNum, width, height);
		createNewToplogicalGraph(newLabels, width, height, newAHGn, objectNum,oNode);

		//展示部分对象信息
		for (int i = 0; i< 10; i++)
			oNode[i].showInformation();

		//19-10-20
		//计算对象序列属性
		//计算vk
		double vk = 0, sumaq = 0, tempaq = 0, spectualStandDeviation = 0;

		for (int i = 0; i<objectNum; i++)
		{
			double tempBsqr = 0, tempGsqr = 0, tempRsqr = 0;
			double tempSqr = 0, sumSqr = 0;
			/*光谱标准差以蓝绿红三个取平均为准*/
			for (int j = 0; j<oNode[i].pixelLocation.size(); j++)
			{
				tempBsqr = (srimg.data[oNode[i].pixelLocation[j]*4] - oNode[i].avgB)*(srimg.data[oNode[i].pixelLocation[j]*4] - oNode[i].avgB);
				tempGsqr = (srimg.data[oNode[i].pixelLocation[j]*4+1] - oNode[i].avgG)*(srimg.data[oNode[i].pixelLocation[j]*4+1] - oNode[i].avgG);
				tempRsqr = (srimg.data[oNode[i].pixelLocation[j]*4+2] - oNode[i].avgR)*(srimg.data[oNode[i].pixelLocation[j]*4+2] - oNode[i].avgR);
				tempSqr = (tempBsqr + tempGsqr + tempRsqr)/3;
				sumSqr += tempSqr;
			}
			spectualStandDeviation = sqrt(sumSqr/oNode[i].pixelnum);
			tempaq = oNode[i].pixelnum*(spectualStandDeviation); //面积乘以光谱标准差
			sumaq += tempaq;
		}
		vk = sumaq/(width*height);
		printf("\n *******************************************\n");
		printf("vk = %lf", vk);
		printf("\n *******************************************\n");

		//计算MI
		double x_meancolor = 0;
		double sumB = 0, sumG = 0, sumR = 0;
		for (int i = 0;i<size.height;i++)
			for (int j = 0; j<size.width; j++)
			{
				sumB += srimg.data[(i*size.width+j)*4];
				sumG += srimg.data[(i*size.width+j)*4+1];
				sumR += srimg.data[(i*size.width+j)*4+2];
			}
		x_meancolor = (sumB + sumG +sumR)/(3*(size.height*size.width));

		int N = 0, sumWij = 0;
		N = objectNum;

		
		for(int i = 0; i<objectNum; i++)
			oNode[i].spectralFeatureInit();  //初始化光谱信息

		double sumXijmean = 0, sumXii = 0;

		for (int i = 0; i<objectNum; i++)
		{
			forward_list<GraphNode>::iterator it;
			for (it = newAHGn[i].pGraphNodeList.begin(); it != newAHGn[i].pGraphNodeList.end(); it++)
			{
				sumWij++; //邻接计数
				sumXijmean += (oNode[i].brightnessBGR - x_meancolor) * (oNode[it->ID].brightnessBGR - x_meancolor);
			}
		}
		
		for (int i = 0; i<objectNum; i++)
		{
			sumXii += (oNode[i].brightnessBGR - x_meancolor)*(oNode[i].brightnessBGR - x_meancolor);
		}

		double MI = 0;
		MI = ((N*sumXijmean)/2) / (sumWij*sumXii);
		printf("\n *******************************************\n");
		printf("MI = %lf", MI);
		printf("\n *******************************************\n");

		//融合效果展示
		Mat imgMerge = srimg.clone();
		for (int i = 1; i<height-1; i++)
			for (int j = 1;j<width-1; j++)
			{
				//不考虑图像边缘
				if (newLabels[i*width + j] != newLabels[(i-1)*width +j] || newLabels[i*width + j] != newLabels[(i+1)*width +j] || newLabels[i*width + j] != newLabels[i*width +j+1] || newLabels[i*width + j] != newLabels[i*width + j-1])
				{
					if (oNode[newLabels[i*width + j]].haveInit == 0)
					{
						oNode[newLabels[i*width + j]].formFeatureInit(width);
						oNode[newLabels[i*width + j]].spectralFeatureInit();
						oNode[newLabels[i*width + j]].haveInit = 1;
					}
					if (oNode[newLabels[i*width+j]].objectTypes == 1)
					{
						imgMerge.data[(i*width+j)*4] = 255;
						imgMerge.data[(i*width+j)*4+1] = 0;
						imgMerge.data[(i*width+j)*4+2] = 0;
						imgMerge.data[(i*width+j)*4+3] = 0;
					}
					else
					{
						imgMerge.data[(i*width+j)*4] = 0;
						imgMerge.data[(i*width+j)*4+1] = 0;
						imgMerge.data[(i*width+j)*4+2] = 255;
						imgMerge.data[(i*width+j)*4+3] = 0;
					}
				}	
			}

			imwrite("out.jpg",imgMerge);
			namedWindow("Superpixel");
			imshow("Superpixel", imgMerge);
			waitKey(0);


	}while(1);


	//超像素分割效果
	Mat imgSuperpixel = srimg.clone();
	for (int i = 1; i<height-1; i++)
		for (int j = 1;j<width-1; j++)
		{
			//不考虑图像边缘
			if (clabels[i*width + j] != clabels[(i-1)*width +j] || clabels[i*width + j] != clabels[(i+1)*width +j] || clabels[i*width + j] != clabels[i*width +j+1] || clabels[i*width + j] != clabels[i*width + j-1])
				{
					imgSuperpixel.data[(i*width+j)*4] = 0;
					imgSuperpixel.data[(i*width+j)*4+1] = 0;
					imgSuperpixel.data[(i*width+j)*4+2] = 255;
					imgSuperpixel.data[(i*width+j)*4+3] = 0;
				}	
		}

	namedWindow("Superpixel");
	imshow("Superpixel", imgSuperpixel);
	waitKey(0);
	//cvCopy(srimg,Imgfinal,NULL);



    //Deallocate memory
    delete rin;
    delete gin;
    delete bin;
    delete lvec;
    delete avec;
    delete bvec;
    delete klabels;
    delete clabels;
    delete seedIndices;
    delete kseedsx;
    delete kseedsy;
    delete kseedsl;
    delete kseedsa;
    delete kseedsb;


	return 0;
}