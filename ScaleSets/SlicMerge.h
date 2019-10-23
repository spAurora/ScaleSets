/*�����ں�*/
/*����������ʵ���ඨ��*/
extern int width;
extern int height;

using namespace std;
class CSuperPixelSet  //������ʵ����
{
public:
	int level;
	int id;
	double avgB;   //BGR��ֵ
	double avgG;
	double avgR;
	double avgNIR;
	int pixelnum;

	vector<int> pixelLocation;  //���ɹ��ɳ����ص�����


	CSuperPixelSet()  //�޲ι���
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

//��������
class BTreeNode
{
public:
	BTreeNode* left;
	BTreeNode* right;

	//������
	int ID;
	//�ȼ���Ϣ
	int level;   //1Ϊ����������
	//������Ϣ
	double avgB;   
	double avgG;
	double avgR;
	double avgNIR;
	//������Ϣ
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
	int ID; //*��ʼ��
	GraphNode(int mID)
	{
		ID = mID;
	}
};


class ArrayHeadGraphNode  //ͷ�������
{
public:
	BTreeNode* pBTnode;							//ָ�������Ľ��
	forward_list<GraphNode> pGraphNodeList;		//�ڽ����˵�
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
	int objectTypes;  //�������
	int haveInit;

	//Э�������
	//��̬ѧ���
	double varxx;
	double varyy;
	double covxy;
	int borderLength;
	double shapeIndex;
	double density;

	//�������
	double brightnessBGRNIR;
	double brightnessBGR;
	double NDVI;
	double NDWI;
	double SBI;
	double BAI;

	ObjectNode()
	{
		borderLength = 0;
		objectTypes = 0;//0δ���࣬1ˮ�壬2ֲ����3������4������5��·
		haveInit = 0;
	}

	//��Ա����

	//չʾ��Ϣ
	void showInformation()
	{
		this->formFeatureInit(width);
		this->spectralFeatureInit();
		printf("ID��%d\n", this->id);
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

	//��ʼ����̬��Ϣ
	//��������Э�������
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

	//��ʼ��������Ϣ
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

//������ĳ���ȼ���ȫ�����
//��ͷ��㿪ʼ���κ�С�ڵ��ڸõȼ��Ľ���Ϊ���õȼ���㡱
void searchTreeNodeWithLevel(BTreeNode* hierarchicalTree, int level, int superPixelNum)
{
	if(hierarchicalTree->level <= level)
		printf("%d -> ",hierarchicalTree->ID);
	if(hierarchicalTree->level > level && hierarchicalTree->left != NULL)
		searchTreeNodeWithLevel(hierarchicalTree->left, level, superPixelNum);
	if(hierarchicalTree->level > level && hierarchicalTree->right != NULL)
		searchTreeNodeWithLevel(hierarchicalTree->right, level, superPixelNum);
};

//����ĳ������л������������ֵ
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

//����ĳ�ȼ��µ����н��
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

//�����¶����㼯��
void createNewObjectSet(int* newLabels, cv::Mat &srimg, ObjectNode* oNode, int objectNum, int width, int height)
{
	for (int i = 0;i<height;i++)
		for (int j = 0;j<width;j++)
		{
			oNode[newLabels[i*width+j]].pixelLocation.push_back(i*width+j);  //�������������������صĹ�ϵ
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
			oNode[i].objectTypes = 0; //���������Ϊδ����
			//oNode[i].formFeatureInit(width);
		}
}

//�����µ��ڽ�ͼ
//ͬʱ�����������߳�
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
				int check = 0; //0Ϊ������
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
				int check = 0; //0Ϊ������
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




//����������ʵ�弯���Լ�������ײ���
//������ͼ�㡢���ߡ������ؽ���顢ԭͼ�񡢳������������������
void createSuperPixelVector(int* label, int width, int height, CSuperPixelSet *csps, cv::Mat& srimg, int superPixelNum, BTreeNode* hierarchicalTree)
{
		
	//����

	for (int i = 0;i<height;i++)
		for (int j = 0;j<width;j++)
		{
			csps[label[i*width+j]].pixelLocation.push_back(i*width+j);  //�������������������صĹ�ϵ
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

//������ʼ������ͼ
//������ͼ�㡢ͼƬ���ߡ�ͼͷ������顢��������Ŀ
void createToplogicalGraph(int*clabels, int width, int height, ArrayHeadGraphNode* mAhgn ,int superPixelnum)
{

	for (int i = 0; i<height - 1; i++)
		for (int j = 0; j < width - 1; j++)
		{
			if (clabels[i*width + j] != clabels[i*width + j + 1])
			{
				forward_list<GraphNode>::iterator it;
				int check = 0; //0Ϊ������
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
				int check = 0; //0Ϊ������
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

		//����
		clock_t startTime,endTime; 
		startTime = clock();
		for (int i = 0; i<superPixelnum;i++)
			mAhgn[i].pGraphNodeList.sort(cmp);
		endTime = clock();
		cout << "����Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
		
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

//������ײ�
double calculateDifference(BTreeNode t1, BTreeNode t2)
{
	return(sqrt( (t1.avgB-t2.avgB)*(t1.avgB-t2.avgB) + (t1.avgG-t2.avgG)*(t1.avgG-t2.avgG) + (t1.avgR-t2.avgR)*(t1.avgR-t2.avgR) + (t1.avgNIR-t2.avgNIR)));
}

//����Ƿ��ظ�
bool whetherThisValueInTheOtherSet(GraphNode mNode, forward_list<GraphNode> list)
{
	forward_list<GraphNode>::iterator it;
	long startTime, endTime;
	startTime = clock();
	for (it = list.begin(); it != list.end(); it++)
		if (it->ID == mNode.ID)
		{
			//endTime = clock();
			//printf("����ظ�ռ�õ�ms:%ld\n", startTime-endTime);
			return true;
		}
	//endTime = clock();
	//printf("����ظ�ռ�õ�ms:%ld\n", startTime-endTime);
	return false;
}

//ɾ�����
void delNode(ArrayHeadGraphNode* mAhgn, int delInId, int delId_1, int delId_2, int graphAndTreeEnd)
{
	//printf("�ѽ���\n");
	if(mAhgn[delInId].pGraphNodeList.empty())
		printf("�ڽ�Ϊ�գ������߼�����!");

	int sum = 0;
	forward_list<GraphNode>::iterator prev=mAhgn[delInId].pGraphNodeList.before_begin();  //��ʾflst�ġ���ǰԪ�ء�
	forward_list<GraphNode>::iterator curr=mAhgn[delInId].pGraphNodeList.begin();  //��ʾflst�еĵ�һ��Ԫ��
	while(curr!=mAhgn[delInId].pGraphNodeList.end())
	{
		if(curr->ID == delId_1 || curr->ID == delId_2)
		{
			curr=mAhgn[delInId].pGraphNodeList.erase_after(prev);// ɾ�������ƶ�curr
			//printf("��ɾ��\n");
		}
		else
		{
			prev=curr;  //�ƶ�������curr��ָ����һ��Ԫ�أ�prevָ��curr֮ǰ��Ԫ��
			++curr;
		}
	}
	mAhgn[delInId].pGraphNodeList.push_front(GraphNode(graphAndTreeEnd));
}

//���㲢��
void calculateUnion(int childNodeLoc_1, int childNodeLoc_2, int graphAndTreeEnd, ArrayHeadGraphNode* mAhgn, BTreeNode* hierarchicalTree)
{
	mAhgn[childNodeLoc_1].hadRemove = true;
	mAhgn[childNodeLoc_2].hadRemove = true;

	long startTime, endTime;
	startTime = clock();

	//�Ż����ȡ����
	forward_list<GraphNode>::iterator it_1 = mAhgn[childNodeLoc_1].pGraphNodeList.begin();
	forward_list<GraphNode>::iterator it_2 = mAhgn[childNodeLoc_2].pGraphNodeList.begin();
	while (it_1 != mAhgn[childNodeLoc_1].pGraphNodeList.end() && it_2 != mAhgn[childNodeLoc_2].pGraphNodeList.end())
	{
		//Ҫɾ���������ڵ㲻���벢��
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
		else //�������
		{
			mAhgn[graphAndTreeEnd].pGraphNodeList.push_front(*it_2); //������һ��
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
	//��ת
	mAhgn[graphAndTreeEnd].pGraphNodeList.reverse();
	endTime = clock();
	cout << "ȡ����part1:Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

	//check
	if(mAhgn[childNodeLoc_1].pGraphNodeList.empty() == true)
		printf("��1\n");
	if(mAhgn[childNodeLoc_2].pGraphNodeList.empty() == true)
		printf("��2\n");

	//�޸������ڽӽ������˱�
	startTime = clock();
	forward_list<GraphNode>::iterator itt;
	itt = mAhgn[1].pGraphNodeList.begin();
	//printf("WTF %d\n", itt->ID);
	for (itt = mAhgn[childNodeLoc_1].pGraphNodeList.begin(); itt!= mAhgn[childNodeLoc_1].pGraphNodeList.end(); itt++)
	{
		//printf("%d\n", itt->ID);
		if (itt->ID != childNodeLoc_2)  //������˸��ֺţ�������
		{
			//printf("��ʼɾ��2���\n");
			//printf("������ID:it->ID = %d\n", itt->ID);
			delNode(mAhgn, itt->ID, childNodeLoc_1, childNodeLoc_2, graphAndTreeEnd);
		}
	}
	//printf("���1���ڽӵ��������(������Ҫ�ںϵ�������)\n");
	for (itt = mAhgn[childNodeLoc_2].pGraphNodeList.begin(); itt!= mAhgn[childNodeLoc_2].pGraphNodeList.end(); itt++)
	{
		//printf("%d\n", itt->ID);
		if (itt->ID != childNodeLoc_1)
		{
			//printf("��ʼɾ��2���\n");
			//printf("������ID:it->ID = %d\n", itt->ID);
			delNode(mAhgn, itt->ID, childNodeLoc_1, childNodeLoc_2, graphAndTreeEnd);
		}
	}
	endTime = clock();
	cout << "ȡ����part2:Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	//printf("���2���ڽӵ��������(������Ҫ�ںϵ�������)\n");
	/*for (itt = mAhgn[childNodeLoc_2].pGraphNodeList.begin(); itt!= mAhgn[childNodeLoc_2].pGraphNodeList.end(); itt++)
		if (itt->ID != childNodeLoc_1);
		{
			delNode(mAhgn, itt->ID, childNodeLoc_1, childNodeLoc_2, graphAndTreeEnd);
			printf("��ʼɾ��2���\n");
		}*/
}

//�ݹ�����
void DFS(int location,int *vnum, ArrayHeadGraphNode* mAhgn,BTreeNode* hierarchicalTree,int &graphAndTreeEnd, int nowLevel, double & allowDifference, bool& NodeMerge, int superPixelNum, int &mDepth)
{
	if (NodeMerge == true)   //�Ѿ��������ں���ֱ�ӷ���
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
	//------ �ں�
	clock_t startTime,endTime; 
	
	forward_list<GraphNode>::iterator mit;
	for(mit = mAhgn[location].pGraphNodeList.begin(); mit!= mAhgn[location].pGraphNodeList.end(); mit++)
		if (mAhgn[mit->ID].hadRemove == false)
		{
			if(calculateDifference(hierarchicalTree[mit->ID], hierarchicalTree[location]) < allowDifference) //С����ֵ��ں�
			{
				//�ںϰ��� �����½ڵ�
				//ɾ���������ڵ�
				//����Ѿ��������ں�
				//ֱ���˳��ú���
				
				NodeMerge = true;
				graphAndTreeEnd++;
				
				printf("�ںϼ�����ʼ,�����c1 c2�±�Ϊ%d, %d\n, level: %d\n", location, mit->ID, nowLevel);
				//printf("��ʼ�ںϣ�\n");
				startTime = clock();
				calculateUnion(location, mit->ID, graphAndTreeEnd, mAhgn, hierarchicalTree); //����ͼȡ����********
				endTime = clock();
				cout << "ȡ����Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
				//printf("��չ�����±꣺%d\n", graphAndTreeEnd);
				
				//mAhgn[mit->ID].hadRemove = true;
				//mAhgn[location].hadRemove = true;
				//д�������������
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
				//printf("��ֵ����̫�ϸ��޷��ں� code 0\n��ǰ��ֵ��%lf, ������ֵ: %lf\n", calculateDifference(hierarchicalTree[mit->ID], hierarchicalTree[location]), allowDifference);
			}
		}
		else
		{
			//printf("����ѱ��Ƴ����޷��ںϣ�code 1\n");
		}
	
	//------
	//printf("\n��ǰ������ܱ������ڽӵĺ���û���ںϣ�������һ��(�ݹ�):\n**************************************\n");
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


	for(itp = mAhgn[location].pGraphNodeList.begin(); itp != mAhgn[location].pGraphNodeList.end(); itp++)    //mAhgn�����������ᷢ���仯������ʹ�õ�����ʧЧ����������������
	{
		
		//printf("������û��\n");
		//printf("it ID:%d,  vnum[it->ID]:%d,   mAhgn.hadRemove:%d\n", itp->ID, vnum[itp->ID], mAhgn[itp->ID].hadRemove);
		if (vnum[itp->ID] == 0 && mAhgn[itp->ID].hadRemove == false)
		{	int temp;
			temp = itp->ID;
			//printf("temp(Ҳ����it->ID):%d\n");
			DFS(temp, vnum, mAhgn, hierarchicalTree, graphAndTreeEnd, nowLevel, allowDifference, NodeMerge, superPixelNum,mDepth);
		}
		if (NodeMerge == true)
			return;   //!!!!!!!!!!!
	}
}

//�����Լ��ں�
//ͼͷ������顢��������飬ͼ�����������յ㣬��ǰ�����ȼ�������Ĺ��ײ�,�����Ƿ����ں�,��������Ŀ
void traversalAndMerge(ArrayHeadGraphNode* mAhgn,BTreeNode* hierarchicalTree,int &graphAndTreeEnd, int nowLevel, double& allowDifference, bool& NodeMerge, int superPixelNum)  //double�����int��
{
	if (NodeMerge == true)
	{
		printf("something wrong!\n");
		exit(-1);
	}
	int *vnum = new int[2*superPixelNum-1] (); //vnum,0Ϊ��ǰδ������
	for (int i = 0; i<=graphAndTreeEnd;i++)
		if(vnum[i] == 0 && mAhgn[i].hadRemove == false)
		{
			int mDepth = 0;
			//printf("�ӵ�%d��ʼ����\n", i);
			DFS(i, vnum, mAhgn, hierarchicalTree, graphAndTreeEnd, nowLevel, allowDifference, NodeMerge, superPixelNum, mDepth); //ͼ����ȫ��ͨ��
			printf("�ݹ����:%d\n",mDepth);
			break;
		}
	delete vnum;
}

//����ĳ������л������������ֵ
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

//���������
//ͼͷ������顢��������顢ԭͼ�񡢷ֲ�������������
int createHierarchicalTree(ArrayHeadGraphNode* mAhgn,BTreeNode* hierarchicalTree, cv::Mat& srimg, double maxDifference,int superPixelnum,CSuperPixelSet* csps)
{
	int levelindex = 510;  //*��ʱ* ���ڿ��Ʋ���
	int nowLevel = 1;
	bool NodeMerge = false;
	double allowDifference = 0; //����������Բ���
	int graphAndTreeEnd = superPixelnum - 1;
	bool overDifferenceLimit = false;

	FILE *fp;
	if((fp = fopen("Info.txt", "w+")) == NULL)
	{
		printf("��Info��¼�ļ�ʧ��\n");
		exit(-1);
	}

	for (int i = 1; i <= levelindex; i++)
	{
		nowLevel++;

		if (overDifferenceLimit != true) // �����޵�����µ�������
			allowDifference =allowDifference + ((double)510/(double)levelindex); //��ǰΪ�̶�����
		else
			break;

		if(allowDifference >= maxDifference) //����
		{
			allowDifference = 1024;   //���������ֱ�ӽ��������������
			overDifferenceLimit = true;
		}

		printf("allDifference:%lf\n", allowDifference);

		NodeMerge = false;
		do 
		{	
			NodeMerge = false;
			//���������������
			//printf("allDifference:%lf\n", allowDifference);
			//system("pause");
			clock_t startTime,endTime; 
			startTime = clock();
			/*************/
			//�������ںϽڵ㣨ֻ�ں�һ�Σ�
			traversalAndMerge(mAhgn, hierarchicalTree, graphAndTreeEnd, nowLevel, allowDifference, NodeMerge, superPixelnum);
			/*************/
			endTime = clock();
			cout << "Totle Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
			if (NodeMerge == false)
				printf("��ɵ�ǰ��������ֵ�������ں�...\n\n");
		} while (NodeMerge == true);

		//*�ں���ɺ���㵱ǰ��ε�vk��MI��������һ������vk��MI���������Բ���*//

		int level;
		int *newLabels = new int[height*width] ();
		int setvalue = -1;

		for(int i = 0; i <= graphAndTreeEnd; i++) //set all node value
		{
			if (mAhgn[i].hadRemove != true)  // ����ýڵ㲢�����ں�
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
			/*���ױ�׼�������̺�����ȡƽ��Ϊ׼*/
			for (int j = 0; j<oNode[i].pixelLocation.size(); j++)
			{
				tempBsqr = (srimg.data[oNode[i].pixelLocation[j]*4] - oNode[i].avgB)*(srimg.data[oNode[i].pixelLocation[j]*4] - oNode[i].avgB);
				tempGsqr = (srimg.data[oNode[i].pixelLocation[j]*4+1] - oNode[i].avgG)*(srimg.data[oNode[i].pixelLocation[j]*4+1] - oNode[i].avgG);
				tempRsqr = (srimg.data[oNode[i].pixelLocation[j]*4+2] - oNode[i].avgR)*(srimg.data[oNode[i].pixelLocation[j]*4+2] - oNode[i].avgR);
				tempSqr = (tempBsqr + tempGsqr + tempRsqr)/3;
				sumSqr += tempSqr;
			}
			spectualStandDeviation = sqrt(sumSqr/oNode[i].pixelnum);
			tempaq = oNode[i].pixelnum*(spectualStandDeviation); //������Թ��ױ�׼��
			sumaq += tempaq;
		}
		vk = sumaq/(width*height);
		printf("\n *******************************************\n");
		printf("vk = %lf", vk);
		printf("\n *******************************************\n");

		//����MI
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
				oNode[i].spectralFeatureInit();  //��ʼ��������Ϣ

			double sumXijmean = 0, sumXii = 0;

			for (int i = 0; i<objectNum; i++)
			{
				forward_list<GraphNode>::iterator it;
				for (it = newAHGn[i].pGraphNodeList.begin(); it != newAHGn[i].pGraphNodeList.end(); it++)
				{
					sumWij++; //�ڽӼ���
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

		//�ͷ��ڴ�
		free(newLabels);
		printf("�ͷ�newLabels�ɹ�\n");
		//free(oNode);
		//printf("�ͷ�oNode�ɹ�\n");
		//free(newAHGn);
		//printf("�ͷ�newAHGn�ɹ�\n");
	}
	printf("***************\n nowlevel = %d\n*****************\n", nowLevel);
	fclose(fp);
	return nowLevel;
}