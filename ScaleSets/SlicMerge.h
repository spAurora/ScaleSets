/*对象融合*/
/*包括超像素实体类定义*/
extern int width;
extern int height;

using namespace std;
class CSuperPixelSet  //超像素实体类
{
public:
	int level;
	int id;
	double avgB;   //BGR均值
	double avgG;
	double avgR;
	double avgNIR;
	int pixelnum;

	vector<int> pixelLocation;  //容纳构成超像素的像素


	CSuperPixelSet()  //无参构造
	{
		level = 0;
		id = -1;
		avgB = 0;
		avgG = 0;
		avgR = 0;
		avgNIR = 0;
		pixelnum = 0;
	}
	CSuperPixelSet(int flevel, int fid)
	{
		level = flevel;
		id = fid;
	}
protected:
private:
};

//层次树结点
class BTreeNode
{
public:
	BTreeNode* left;
	BTreeNode* right;

	//数据域
	int ID;
	//等级信息
	int level;   //1为基础超像素
	//光谱信息
	double avgB;   
	double avgG;
	double avgR;
	double avgNIR;
	//像素信息
	int pixelnum;

	BTreeNode()
	{
		ID = 0;
		level = 1;
		avgB = 0;
		avgG = 0;
		avgR = 0;
		avgNIR = 0;
		left = NULL;
		right = NULL;
	}

protected:
private:
};

class GraphNode
{
public:
	int ID; //*初始化
	GraphNode(int mID)
	{
		ID = mID;
	}
};


class ArrayHeadGraphNode  //头结点数组
{
public:
	BTreeNode* pBTnode;							//指向层次树的结点
	forward_list<GraphNode> pGraphNodeList;		//邻接拓扑点
	bool hadRemove;

	ArrayHeadGraphNode()
	{
		pBTnode = NULL;
		hadRemove = false; 
	}
protected:
private:
};

class ObjectNode : public CSuperPixelSet
{
public:
	int objectTypes;  //地物类别
	int haveInit;

	//协方差矩阵
	//形态学相关
	double varxx;
	double varyy;
	double covxy;
	int borderLength;
	double shapeIndex;
	double density;

	//光谱相关
	double brightnessBGRNIR;
	double brightnessBGR;
	double NDVI;
	double NDWI;
	double SBI;
	double BAI;

	ObjectNode()
	{
		borderLength = 0;
		objectTypes = 0;//0未分类，1水体，2植被，3裸土，4建筑，5道路
		haveInit = 0;
	}

	//成员函数

	//展示信息
	void showInformation()
	{
		this->formFeatureInit(width);
		this->spectralFeatureInit();
		printf("ID：%d\n", this->id);
		printf("Area: %d\n", this->pixelnum);
		printf("borderLength:%d\n", this->borderLength);
		printf("shapeIndex: %lf\n", this->shapeIndex);
		printf("density:%lf\n", this->density);
		printf("brightnessBGRNIR:%lf\n", this->brightnessBGRNIR);
		printf("brightnessBGR:%lf\n", this->brightnessBGR);
		printf("NDVI:%lf\n", this->NDVI);
		printf("NDWI:%lf\n", this->NDWI);
		printf("BAI:%lf\n", this->BAI);
		printf("SBI:%lf\n\n", this->SBI);

	}

	//初始化形态信息
	//计算依赖协方差矩阵
	void formFeatureInit(int width)
	{
		vector<int>::iterator it;
		int sumx = 0;
		int sumy = 0;
		for(it = this->pixelLocation.begin(); it != this->pixelLocation.end(); it++)
		{
			sumx += (*it)/width;
			sumy += (*it)%width;
		}
		double meanx = 0, meany = 0;
		meanx = sumx/(double)this->pixelLocation.size();
		meany = sumy/(double)this->pixelLocation.size();
		double finx = 0, finy = 0, finxy = 0;
		for(it = this->pixelLocation.begin(); it != this->pixelLocation.end(); it++)
		{
			finx += ((*it)/width - meanx)*((*it)/width - meanx);
			finy += ((*it)%width - meany)*((*it)%width - meany);
			finxy += ((*it)/width - meanx) * ((*it)%width - meany);
		}
		this->varxx = finx/(this->pixelLocation.size()-1);
		this->varyy = finy/(this->pixelLocation.size()-1);
		this->covxy = finxy/(this->pixelLocation.size()-1);

		this->shapeIndex = (double)borderLength/(4*sqrt((double)this->pixelnum));
		this->density = sqrt((double)this->pixelnum)/(1+sqrt(varxx+varyy));
	}

	//初始化光谱信息
	void spectralFeatureInit()
	{
		this->NDVI = (this->avgNIR - this->avgR)/(double)(this->avgNIR + this->avgR);
		this->NDWI = (this->avgG - this->avgNIR)/(double)(this->avgG + this->avgNIR);
		this->brightnessBGRNIR = (this->avgB + this->avgG + this->avgR + this->avgNIR)/4;
		this->brightnessBGR = (this->avgB + this->avgG + this->avgR)/3;
		this->BAI = (this->avgB - this ->avgNIR)/(this->avgB + this->avgNIR);
		this->SBI = sqrt(this->avgR*this->avgR + this->avgNIR*this->avgNIR);
	}



protected:
private:
};

//遍历出某个等级的全部结点
//从头结点开始，任何小于等于该等级的结点均为“该等级结点”
void searchTreeNodeWithLevel(BTreeNode* hierarchicalTree, int level, int superPixelNum)
{
	if(hierarchicalTree->level <= level)
		printf("%d -> ",hierarchicalTree->ID);
	if(hierarchicalTree->level > level && hierarchicalTree->left != NULL)
		searchTreeNodeWithLevel(hierarchicalTree->left, level, superPixelNum);
	if(hierarchicalTree->level > level && hierarchicalTree->right != NULL)
		searchTreeNodeWithLevel(hierarchicalTree->right, level, superPixelNum);
};

//设置某结点所有基层结点的子像素值
void setNowLevelNodeValue(int *labels, int &setValue, BTreeNode* htNode, CSuperPixelSet* csps)
{
	if (htNode->level == 1)
	{
		vector<int>::iterator it;
		for(it = csps[htNode->ID].pixelLocation.begin(); it != csps[htNode->ID].pixelLocation.end(); it++)
			labels[*it] = setValue;
	}
	if(htNode->level > 1 && htNode->left != NULL)
		setNowLevelNodeValue(labels, setValue, htNode->left, csps);
	if(htNode->level > 1 && htNode->right != NULL)
		setNowLevelNodeValue(labels, setValue, htNode->right, csps);
}

//遍历某等级下的所有结点
void setAllNodeValue(int *labels, int level, BTreeNode* hierarchicalTree, int &setValue, CSuperPixelSet* csps)
{
	if(hierarchicalTree->level <= level)
	{
		setValue++;
		setNowLevelNodeValue(labels, setValue, hierarchicalTree, csps);
	}
	if(hierarchicalTree->level > level && hierarchicalTree->left != NULL)
		setAllNodeValue(labels, level, hierarchicalTree->left, setValue, csps);
	if(hierarchicalTree->level > level && hierarchicalTree->right != NULL)
		setAllNodeValue(labels, level, hierarchicalTree->right, setValue, csps);
}

//建立新对象结点集合
void createNewObjectSet(int* newLabels, cv::Mat &srimg, ObjectNode* oNode, int objectNum, int width, int height)
{
	for (int i = 0;i<height;i++)
		for (int j = 0;j<width;j++)
		{
			oNode[newLabels[i*width+j]].pixelLocation.push_back(i*width+j);  //建立超像素与其中像素的关系
			oNode[newLabels[i*width+j]].pixelnum++;
			oNode[newLabels[i*width+j]].avgB += srimg.data[(i*width+j)*4];
			oNode[newLabels[i*width+j]].avgG += srimg.data[(i*width+j)*4+1];
			oNode[newLabels[i*width+j]].avgR += srimg.data[(i*width+j)*4+2];
			oNode[newLabels[i*width+j]].avgNIR += srimg.data[(i*width+j)*4+3];
		}
		for(int i = 0; i< objectNum;i++)
		{
			oNode[i].avgB /= oNode[i].pixelnum;
			oNode[i].avgG /= oNode[i].pixelnum;
			oNode[i].avgR /= oNode[i].pixelnum;
			oNode[i].avgNIR /= oNode[i].pixelnum;
			oNode[i].id = i;
			oNode[i].objectTypes = 0; //定义地物类为未定义
			//oNode[i].formFeatureInit(width);
		}
}

//建立新的邻接图
//同时计算各个对象边长
void createNewToplogicalGraph(int *newLabels, int width, int height, ArrayHeadGraphNode* newAHGn, int objectNum, ObjectNode* oNode)
{
	for (int i = 0; i<height; i++)
		oNode[newLabels[i*width+width-1]].borderLength++;
	for(int j = 0; j<width-1; j++)
		oNode[newLabels[(height-1)*width+j]].borderLength++;
	for (int i = 0; i<height - 1; i++)
		for (int j = 0; j < width - 1; j++)
		{
			if(i == 0 || j == 0)
				oNode[newLabels[i*width+j]].borderLength++;
			if (newLabels[i*width + j] != newLabels[i*width + j + 1])
			{
				oNode[newLabels[i*width+j]].borderLength++;
				forward_list<GraphNode>::iterator it;
				int check = 0; //0为不存在
				for (it = newAHGn[newLabels[i*width + j]].pGraphNodeList.begin(); it != newAHGn[newLabels[i*width + j]].pGraphNodeList.end(); it++)
					if (it->ID == newLabels[i*width+j+1])
					{
						check = 1;
						break;
					}
					if(check == 0)
					{
						newAHGn[newLabels[i*width + j]].pGraphNodeList.push_front(newLabels[i*width + j + 1]);
						newAHGn[newLabels[i*width + j + 1]].pGraphNodeList.push_front(newLabels[i*width + j]);
					}
			}

			if (newLabels[i*width + j] != newLabels[(i+1)*width+j])
			{
				oNode[newLabels[i*width+j]].borderLength++;
				forward_list<GraphNode>::iterator it;
				int check = 0; //0为不存在
				for (it = newAHGn[newLabels[i*width + j]].pGraphNodeList.begin(); it != newAHGn[newLabels[i*width + j]].pGraphNodeList.end(); it++)
					if (it->ID == newLabels[(i+1)*width+j])
					{
						check = 1;
						break;
					}
					if(check == 0)
					{
						newAHGn[newLabels[i*width + j]].pGraphNodeList.push_front(newLabels[(i+1)*width+j]);
						newAHGn[newLabels[(i+1)*width+j]].pGraphNodeList.push_front(newLabels[i*width + j]);
					}
			}
		}
}




//创建超像素实体集合以及层次树底层结点
//超像素图层、宽、高、超像素结点组、原图像、超像素数、层次数数组
void createSuperPixelVector(int* label, int width, int height, CSuperPixelSet *csps, cv::Mat& srimg, int superPixelNum, BTreeNode* hierarchicalTree)
{
		
	//排序

	for (int i = 0;i<height;i++)
		for (int j = 0;j<width;j++)
		{
			csps[label[i*width+j]].pixelLocation.push_back(i*width+j);  //建立超像素与其中像素的关系
			csps[label[i*width+j]].pixelnum++;
			csps[label[i*width+j]].avgB += srimg.data[(i*width+j)*4];
			csps[label[i*width+j]].avgG += srimg.data[(i*width+j)*4+1];
			csps[label[i*width+j]].avgR += srimg.data[(i*width+j)*4+2];
			csps[label[i*width+j]].avgNIR += srimg.data[(i*width+j)*4+3];
		}
	for(int i = 0; i< superPixelNum;i++)
	{
		csps[i].avgB /= csps[i].pixelnum;
		csps[i].avgG /= csps[i].pixelnum;
		csps[i].avgR /= csps[i].pixelnum;
		csps[i].avgNIR /= csps[i].pixelnum;
	}
	printf("-----------\n %lf \n ----------", csps[100].avgB);

	for (int i = 0; i< superPixelNum; i++)
	{
		hierarchicalTree[i].level = 1;
		hierarchicalTree[i].ID = i;
		hierarchicalTree[i].avgB = csps[i].avgB;
		hierarchicalTree[i].avgG = csps[i].avgG;
		hierarchicalTree[i].avgR = csps[i].avgR;
		hierarchicalTree[i].avgNIR = csps[i].avgNIR;
		hierarchicalTree[i].pixelnum = csps[i].pixelnum;
	}
}

bool cmp(GraphNode first, GraphNode second) {
	return first.ID > second.ID;
}

//创建初始化拓扑图
//超像素图层、图片宽、高、图头结点数组、超像素数目
void createToplogicalGraph(int*clabels, int width, int height, ArrayHeadGraphNode* mAhgn ,int superPixelnum)
{

	for (int i = 0; i<height - 1; i++)
		for (int j = 0; j < width - 1; j++)
		{
			if (clabels[i*width + j] != clabels[i*width + j + 1])
			{
				forward_list<GraphNode>::iterator it;
				int check = 0; //0为不存在
				for (it = mAhgn[clabels[i*width + j]].pGraphNodeList.begin(); it != mAhgn[clabels[i*width + j]].pGraphNodeList.end(); it++)
					if (it->ID == clabels[i*width+j+1])
					{
						check = 1;
						break;
					}
				if(check == 0)
				{
					mAhgn[clabels[i*width + j]].pGraphNodeList.push_front(clabels[i*width + j + 1]);
					mAhgn[clabels[i*width + j + 1]].pGraphNodeList.push_front(clabels[i*width + j]);
				}
			}

			if (clabels[i*width + j] != clabels[(i+1)*width+j])
			{
				forward_list<GraphNode>::iterator it;
				int check = 0; //0为不存在
				for (it = mAhgn[clabels[i*width + j]].pGraphNodeList.begin(); it != mAhgn[clabels[i*width + j]].pGraphNodeList.end(); it++)
					if (it->ID == clabels[(i+1)*width+j])
					{
						check = 1;
						break;
					}
					if(check == 0)
					{
						mAhgn[clabels[i*width + j]].pGraphNodeList.push_front(clabels[(i+1)*width+j]);
						mAhgn[clabels[(i+1)*width+j]].pGraphNodeList.push_front(clabels[i*width + j]);
					}
			}
		}

		//排序
		clock_t startTime,endTime; 
		startTime = clock();
		for (int i = 0; i<superPixelnum;i++)
			mAhgn[i].pGraphNodeList.sort(cmp);
		endTime = clock();
		cout << "排序Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
		
		forward_list<GraphNode>::iterator it;
		for (it = mAhgn[30].pGraphNodeList.begin(); it!= mAhgn[30].pGraphNodeList.end(); it++)
			printf("%d -> ", it->ID);
		printf("\n");
		for (it = mAhgn[256].pGraphNodeList.begin(); it!= mAhgn[256].pGraphNodeList.end(); it++)
			printf("%d -> ", it->ID);
		printf("\n");
		for (it = mAhgn[1024].pGraphNodeList.begin(); it!= mAhgn[1024].pGraphNodeList.end(); it++)
			printf("%d -> ", it->ID);
		printf("\n");
}

//计算光谱差
double calculateDifference(BTreeNode t1, BTreeNode t2)
{
	return(sqrt( (t1.avgB-t2.avgB)*(t1.avgB-t2.avgB) + (t1.avgG-t2.avgG)*(t1.avgG-t2.avgG) + (t1.avgR-t2.avgR)*(t1.avgR-t2.avgR) + (t1.avgNIR-t2.avgNIR)));
}

//检查是否重复
bool whetherThisValueInTheOtherSet(GraphNode mNode, forward_list<GraphNode> list)
{
	forward_list<GraphNode>::iterator it;
	long startTime, endTime;
	startTime = clock();
	for (it = list.begin(); it != list.end(); it++)
		if (it->ID == mNode.ID)
		{
			//endTime = clock();
			//printf("检查重复占用的ms:%ld\n", startTime-endTime);
			return true;
		}
	//endTime = clock();
	//printf("检查重复占用的ms:%ld\n", startTime-endTime);
	return false;
}

//删除结点
void delNode(ArrayHeadGraphNode* mAhgn, int delInId, int delId_1, int delId_2, int graphAndTreeEnd)
{
	//printf("已进入\n");
	if(mAhgn[delInId].pGraphNodeList.empty())
		printf("邻接为空，存在逻辑错误!");

	int sum = 0;
	forward_list<GraphNode>::iterator prev=mAhgn[delInId].pGraphNodeList.before_begin();  //表示flst的“首前元素”
	forward_list<GraphNode>::iterator curr=mAhgn[delInId].pGraphNodeList.begin();  //表示flst中的第一个元素
	while(curr!=mAhgn[delInId].pGraphNodeList.end())
	{
		if(curr->ID == delId_1 || curr->ID == delId_2)
		{
			curr=mAhgn[delInId].pGraphNodeList.erase_after(prev);// 删除它并移动curr
			//printf("已删除\n");
		}
		else
		{
			prev=curr;  //移动迭代器curr，指向下一个元素，prev指向curr之前的元素
			++curr;
		}
	}
	mAhgn[delInId].pGraphNodeList.push_front(GraphNode(graphAndTreeEnd));
}

//计算并集
void calculateUnion(int childNodeLoc_1, int childNodeLoc_2, int graphAndTreeEnd, ArrayHeadGraphNode* mAhgn, BTreeNode* hierarchicalTree)
{
	mAhgn[childNodeLoc_1].hadRemove = true;
	mAhgn[childNodeLoc_2].hadRemove = true;

	long startTime, endTime;
	startTime = clock();

	//优化后的取并集
	forward_list<GraphNode>::iterator it_1 = mAhgn[childNodeLoc_1].pGraphNodeList.begin();
	forward_list<GraphNode>::iterator it_2 = mAhgn[childNodeLoc_2].pGraphNodeList.begin();
	while (it_1 != mAhgn[childNodeLoc_1].pGraphNodeList.end() && it_2 != mAhgn[childNodeLoc_2].pGraphNodeList.end())
	{
		//要删除的两个节点不进入并集
		if (it_1->ID == childNodeLoc_1 || it_1->ID == childNodeLoc_2)
		{
			it_1++;
			continue;
		}
		if (it_2->ID == childNodeLoc_1 || it_2->ID == childNodeLoc_2)
		{
			it_2++;
			continue;
		}
		if(it_1->ID > it_2->ID)
		{
			mAhgn[graphAndTreeEnd].pGraphNodeList.push_front(*it_1);
			it_1++;
		}
		else if(it_1->ID < it_2->ID)
		{
			mAhgn[graphAndTreeEnd].pGraphNodeList.push_front(*it_2);
			it_2++;
		}
		else //两者相等
		{
			mAhgn[graphAndTreeEnd].pGraphNodeList.push_front(*it_2); //随便进入一个
			it_1++;
			it_2++;
		}
	}
	while(it_1 != mAhgn[childNodeLoc_1].pGraphNodeList.end())
	{
		mAhgn[graphAndTreeEnd].pGraphNodeList.push_front(*it_1);
		it_1++;
	}
	while(it_2 != mAhgn[childNodeLoc_2].pGraphNodeList.end())
	{
		mAhgn[graphAndTreeEnd].pGraphNodeList.push_front(*it_2);
		it_2++;
	}
	//翻转
	mAhgn[graphAndTreeEnd].pGraphNodeList.reverse();
	endTime = clock();
	cout << "取并集part1:Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

	//check
	if(mAhgn[childNodeLoc_1].pGraphNodeList.empty() == true)
		printf("空1\n");
	if(mAhgn[childNodeLoc_2].pGraphNodeList.empty() == true)
		printf("空2\n");

	//修改其余邻接结点的拓扑表
	startTime = clock();
	forward_list<GraphNode>::iterator itt;
	itt = mAhgn[1].pGraphNodeList.begin();
	//printf("WTF %d\n", itt->ID);
	for (itt = mAhgn[childNodeLoc_1].pGraphNodeList.begin(); itt!= mAhgn[childNodeLoc_1].pGraphNodeList.end(); itt++)
	{
		//printf("%d\n", itt->ID);
		if (itt->ID != childNodeLoc_2)  //后面多了个分号，见鬼了
		{
			//printf("开始删除2结点\n");
			//printf("清理结点ID:it->ID = %d\n", itt->ID);
			delNode(mAhgn, itt->ID, childNodeLoc_1, childNodeLoc_2, graphAndTreeEnd);
		}
	}
	//printf("结点1的邻接点清理完毕(不包括要融合的两个点)\n");
	for (itt = mAhgn[childNodeLoc_2].pGraphNodeList.begin(); itt!= mAhgn[childNodeLoc_2].pGraphNodeList.end(); itt++)
	{
		//printf("%d\n", itt->ID);
		if (itt->ID != childNodeLoc_1)
		{
			//printf("开始删除2结点\n");
			//printf("清理结点ID:it->ID = %d\n", itt->ID);
			delNode(mAhgn, itt->ID, childNodeLoc_1, childNodeLoc_2, graphAndTreeEnd);
		}
	}
	endTime = clock();
	cout << "取并集part2:Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	//printf("结点2的邻接点清理完毕(不包括要融合的两个点)\n");
	/*for (itt = mAhgn[childNodeLoc_2].pGraphNodeList.begin(); itt!= mAhgn[childNodeLoc_2].pGraphNodeList.end(); itt++)
		if (itt->ID != childNodeLoc_1);
		{
			delNode(mAhgn, itt->ID, childNodeLoc_1, childNodeLoc_2, graphAndTreeEnd);
			printf("开始删除2结点\n");
		}*/
}

//递归深搜
void DFS(int location,int *vnum, ArrayHeadGraphNode* mAhgn,BTreeNode* hierarchicalTree,int &graphAndTreeEnd, int nowLevel, double & allowDifference, bool& NodeMerge, int superPixelNum, int &mDepth)
{
	if (NodeMerge == true)   //已经发生过融合则直接返回
		return;
	if (location == 2*superPixelNum-2)
	{
		NodeMerge = false;
		return;
	}
	if(vnum[location]!=0)
		return;
	mDepth++;
	vnum[location] = 1;
	//printf("%d\n", location);
	//------ 融合
	clock_t startTime,endTime; 
	
	forward_list<GraphNode>::iterator mit;
	for(mit = mAhgn[location].pGraphNodeList.begin(); mit!= mAhgn[location].pGraphNodeList.end(); mit++)
		if (mAhgn[mit->ID].hadRemove == false)
		{
			if(calculateDifference(hierarchicalTree[mit->ID], hierarchicalTree[location]) < allowDifference) //小于阈值差，融合
			{
				//融合包括 创建新节点
				//删除这两个节点
				//标记已经发生过融合
				//直接退出该函数
				
				NodeMerge = true;
				graphAndTreeEnd++;
				
				printf("融合即将开始,两结点c1 c2下标为%d, %d\n, level: %d\n", location, mit->ID, nowLevel);
				//printf("开始融合：\n");
				startTime = clock();
				calculateUnion(location, mit->ID, graphAndTreeEnd, mAhgn, hierarchicalTree); //拓扑图取并集********
				endTime = clock();
				cout << "取并集Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
				//printf("扩展到的下标：%d\n", graphAndTreeEnd);
				
				//mAhgn[mit->ID].hadRemove = true;
				//mAhgn[location].hadRemove = true;
				//写入新树结点数据
				hierarchicalTree[graphAndTreeEnd].left = &hierarchicalTree[mit->ID];
				hierarchicalTree[graphAndTreeEnd].right = &hierarchicalTree[location];
				hierarchicalTree[graphAndTreeEnd].pixelnum = hierarchicalTree[location].pixelnum + hierarchicalTree[mit->ID].pixelnum;
				hierarchicalTree[graphAndTreeEnd].level = nowLevel;
				hierarchicalTree[graphAndTreeEnd].avgB = (hierarchicalTree[location].pixelnum*hierarchicalTree[location].avgB + hierarchicalTree[mit->ID].pixelnum*hierarchicalTree[mit->ID].avgB) / hierarchicalTree[graphAndTreeEnd].pixelnum;
				hierarchicalTree[graphAndTreeEnd].avgG = (hierarchicalTree[location].pixelnum*hierarchicalTree[location].avgG + hierarchicalTree[mit->ID].pixelnum*hierarchicalTree[mit->ID].avgG) / hierarchicalTree[graphAndTreeEnd].pixelnum;
				hierarchicalTree[graphAndTreeEnd].avgR = (hierarchicalTree[location].pixelnum*hierarchicalTree[location].avgR + hierarchicalTree[mit->ID].pixelnum*hierarchicalTree[mit->ID].avgR) / hierarchicalTree[graphAndTreeEnd].pixelnum;
				hierarchicalTree[graphAndTreeEnd].avgNIR = (hierarchicalTree[location].pixelnum*hierarchicalTree[location].avgNIR + hierarchicalTree[mit->ID].pixelnum*hierarchicalTree[mit->ID].avgNIR) / hierarchicalTree[graphAndTreeEnd].pixelnum;
				hierarchicalTree[graphAndTreeEnd].ID = graphAndTreeEnd;
				
				return;
			}
			else
			{
				//printf("阈值限制太严格，无法融合 code 0\n当前阈值：%lf, 允许阈值: %lf\n", calculateDifference(hierarchicalTree[mit->ID], hierarchicalTree[location]), allowDifference);
			}
		}
		else
		{
			//printf("结点已被移除，无法融合：code 1\n");
		}
	
	//------
	//printf("\n当前结点与周边所有邻接的好像都没法融合，进行下一步(递归):\n**************************************\n");
	//forward_list<GraphNode>::iterator itt;
	//itt = mAhgn[1].pGraphNodeList.begin();
	////printf("WTF %d\n", itt->ID);
	//for (itt = mAhgn[childNodeLoc_1].pGraphNodeList.begin(); itt!= mAhgn[childNodeLoc_1].pGraphNodeList.end(); itt++)
	forward_list<GraphNode>::iterator itp;
	itp = mAhgn[0].pGraphNodeList.begin();
	//printf("WTF %d\n", itp->ID);
	//for(it = mAhgn[location].pGraphNodeList.begin(); it!= mAhgn[location].pGraphNodeList.end(); it++)
	//for(itp = mAhgn[location].pGraphNodeList.begin(); itp != mAhgn[location].pGraphNodeList.end(); itp++)
	//{
		//printf("FUCK %d\n", itp->ID);
	//}


	for(itp = mAhgn[location].pGraphNodeList.begin(); itp != mAhgn[location].pGraphNodeList.end(); itp++)    //mAhgn的链表容器会发生变化，进而使得迭代器失效！！！！！！！！
	{
		
		//printf("进来了没？\n");
		//printf("it ID:%d,  vnum[it->ID]:%d,   mAhgn.hadRemove:%d\n", itp->ID, vnum[itp->ID], mAhgn[itp->ID].hadRemove);
		if (vnum[itp->ID] == 0 && mAhgn[itp->ID].hadRemove == false)
		{	int temp;
			temp = itp->ID;
			//printf("temp(也就是it->ID):%d\n");
			DFS(temp, vnum, mAhgn, hierarchicalTree, graphAndTreeEnd, nowLevel, allowDifference, NodeMerge, superPixelNum,mDepth);
		}
		if (NodeMerge == true)
			return;   //!!!!!!!!!!!
	}
}

//遍历以及融合
//图头结点数组、层次树数组，图与层次数遍历终点，当前树结点等级，允许的光谱差,遍历是否发生融合,超像素数目
void traversalAndMerge(ArrayHeadGraphNode* mAhgn,BTreeNode* hierarchicalTree,int &graphAndTreeEnd, int nowLevel, double& allowDifference, bool& NodeMerge, int superPixelNum)  //double定义成int了
{
	if (NodeMerge == true)
	{
		printf("something wrong!\n");
		exit(-1);
	}
	int *vnum = new int[2*superPixelNum-1] (); //vnum,0为当前未遍历到
	for (int i = 0; i<=graphAndTreeEnd;i++)
		if(vnum[i] == 0 && mAhgn[i].hadRemove == false)
		{
			int mDepth = 0;
			//printf("从第%d开始遍历\n", i);
			DFS(i, vnum, mAhgn, hierarchicalTree, graphAndTreeEnd, nowLevel, allowDifference, NodeMerge, superPixelNum, mDepth); //图是完全连通的
			printf("递归深度:%d\n",mDepth);
			break;
		}
	delete vnum;
}

//设置某结点所有基层结点的子像素值
void setNowLevelNodeValue_InSlicMerge(int *labels, int &setValue, BTreeNode* htNode, CSuperPixelSet* csps)
{
	if (htNode->level == 1)
	{
		vector<int>::iterator it;
		for(it = csps[htNode->ID].pixelLocation.begin(); it != csps[htNode->ID].pixelLocation.end(); it++)
			labels[*it] = setValue;
	}
	if(htNode->level > 1 && htNode->left != NULL)
		setNowLevelNodeValue_InSlicMerge(labels, setValue, htNode->left, csps);
	if(htNode->level > 1 && htNode->right != NULL)
		setNowLevelNodeValue_InSlicMerge(labels, setValue, htNode->right, csps);
}

//创建层次树
//图头结点数组、层次树数组、原图像、分层数、超像素数
int createHierarchicalTree(ArrayHeadGraphNode* mAhgn,BTreeNode* hierarchicalTree, cv::Mat& srimg, double maxDifference,int superPixelnum,CSuperPixelSet* csps)
{
	int levelindex = 510;  //*临时* 用于控制步长
	int nowLevel = 1;
	bool NodeMerge = false;
	double allowDifference = 0; //允许的异质性差异
	int graphAndTreeEnd = superPixelnum - 1;
	bool overDifferenceLimit = false;

	FILE *fp;
	if((fp = fopen("Info.txt", "w+")) == NULL)
	{
		printf("打开Info记录文件失败\n");
		exit(-1);
	}

	for (int i = 1; i <= levelindex; i++)
	{
		nowLevel++;

		if (overDifferenceLimit != true) // 不超限的情况下调整步长
			allowDifference =allowDifference + ((double)510/(double)levelindex); //当前为固定步长
		else
			break;

		if(allowDifference >= maxDifference) //超限
		{
			allowDifference = 1024;   //超限情况下直接将允许差调整到最大
			overDifferenceLimit = true;
		}

		printf("allDifference:%lf\n", allowDifference);

		NodeMerge = false;
		do 
		{	
			NodeMerge = false;
			//遍历并创建层次树
			//printf("allDifference:%lf\n", allowDifference);
			//system("pause");
			clock_t startTime,endTime; 
			startTime = clock();
			/*************/
			//遍历并融合节点（只融合一次）
			traversalAndMerge(mAhgn, hierarchicalTree, graphAndTreeEnd, nowLevel, allowDifference, NodeMerge, superPixelnum);
			/*************/
			endTime = clock();
			cout << "Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
			if (NodeMerge == false)
				printf("完成当前异质性阈值下所有融合...\n\n");
		} while (NodeMerge == true);

		//*融合完成后计算当前层次的vk和MI，并在下一层依据vk和MI调整异质性步长*//

		int level;
		int *newLabels = new int[height*width] ();
		int setvalue = -1;

		for(int i = 0; i <= graphAndTreeEnd; i++) //set all node value
		{
			if (mAhgn[i].hadRemove != true)  // 如果该节点并参与融合
			{
				setvalue++;
				setNowLevelNodeValue_InSlicMerge(newLabels, setvalue, &hierarchicalTree[i], csps);
			}
		}
		
		int objectNum = setvalue + 1;
		ObjectNode* oNode = new ObjectNode[objectNum];
		ArrayHeadGraphNode *newAHGn = new ArrayHeadGraphNode[objectNum];
		createNewObjectSet(newLabels, srimg, oNode, objectNum, width, height);
		createNewToplogicalGraph(newLabels, width, height, newAHGn, objectNum, oNode);

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
		for (int i = 0;i<height;i++)
			for (int j = 0; j<width; j++)
			{
				sumB += srimg.data[(i*width+j)*4];
				sumG += srimg.data[(i*width+j)*4+1];
				sumR += srimg.data[(i*width+j)*4+2];
			}
			x_meancolor = (sumB + sumG +sumR)/(3*(height*width));

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
		
		fprintf(fp, "%d %d %lf %lf\n", i, objectNum, vk, MI);

		//释放内存
		free(newLabels);
		printf("释放newLabels成功\n");
		//free(oNode);
		//printf("释放oNode成功\n");
		//free(newAHGn);
		//printf("释放newAHGn成功\n");
	}
	printf("***************\n nowlevel = %d\n*****************\n", nowLevel);
	fclose(fp);
	return nowLevel;
}