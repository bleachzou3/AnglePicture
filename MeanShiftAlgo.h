#ifndef MEAN_SHIFT_ALGO_HPP
#define MEAN_SHIFT_ALGO_HPP
#include<map>
#include<list>
#include "MeanPoint.h"
#include <itkImage.h>
#include <stdlib.h>
#include <set>
#include <climits>
#include <log4cpp/Category.hh>
#include <log4cpp/PropertyConfigurator.hh>
#include "CommonUtility.h"

extern log4cpp::Category& rootLog;
extern log4cpp::Category& subLog ;
struct PixelPointCmp
{
	bool operator()(PixelPoint*p1,PixelPoint*p2)
	{
		if(p1->X != p2->X)
		{
			return p1->X > p2->X;
		}
		if(p1->Y != p2->Y)
		{
			return p1->Y > p2->Y;
		}
		return p1->Z > p2->Z;
	}
};
struct MeanPointCmp
{
	bool operator()(MeanPoint*p1,MeanPoint*p2)
	{
		if(p1->X != p2->X)
		{
			return p1->X > p2->X;
		}
		if(p1->Y != p2->Y)
		{
			return p1->Y > p2->Y;
		}
		return p1->Z > p2->Z;
	}
};
typedef set<PixelPoint*,PixelPointCmp> SetPixelPoint;
typedef map<MeanPoint*,SetPixelPoint,MeanPointCmp> MapPixel;
class MeanShiftAlgo
{
private:
	int radius;//点的空间位置之间的距离
	float grayDistance;//点的灰度值之间的距离
	int numOfCenters;//每一轮迭代的时候，求密度最大的点的移动方向
	int iterationNums;//总共需要的迭代次数
	MapPixel allPoints;
	int clusterPointSize;//当合并完成的时候，一个cluster至少的size，如果少于这个size，这个cluster要被归并到其它簇类里面
public:
	MeanShiftAlgo(int _radius,float _grayDistance, int _numOfCenters,int _itertionNums,int _clusterPointSize)
	{
		this->radius = _radius;
		this->grayDistance = _grayDistance;
		this->numOfCenters = _numOfCenters;
		this->iterationNums = _itertionNums;
		this->clusterPointSize = _clusterPointSize;
		
	}
	void filter(itk::Image<signed short,3>*source,itk::Image<signed short,3>*dest)
	{
		subLog.info("void MeanShiftAlgo::filter()开始执行");
	    rootLog.info("void MeanShiftAlgo::filter()开始执行");


		itk::Image<signed short,3>::SizeType size = source->GetLargestPossibleRegion().GetSize();
		itk::Image<signed short,3>::IndexType startIndex = source->GetLargestPossibleRegion().GetIndex();
		itk::Image<signed short,3>::IndexType upperIndex = source->GetLargestPossibleRegion().GetUpperIndex();
		subLog.info("void MeanShiftAlgo::filter()获取图像维度成功"+toString(size[0])+" "+toString(size[1])+" "+toString(size[2])+"startIndex  "+toString(startIndex[0])+" "+toString(startIndex[1])+" "+toString(startIndex[2])+" upperIndex"+toString(upperIndex[0])+" "+toString(upperIndex[1])+" "+toString(upperIndex[2]));
		rootLog.info("void MeanShiftAlgo::filter()获取图像维度成功"+toString(size[0])+" "+toString(size[1])+" "+toString(size[2])+"startIndex  "+toString(startIndex[0])+" "+toString(startIndex[1])+" "+toString(startIndex[2])+" upperIndex"+toString(upperIndex[0])+" "+toString(upperIndex[1])+" "+toString(upperIndex[2]));
		vector<vector<vector<float>>> pixels(size[0],vector<vector<float>>(size[1],vector<float>(size[2],0)));
		vector<vector<vector<bool>>> pixelFlag(size[0],vector<vector<bool>>(size[1],vector<bool>(size[2],false)));
		subLog.info("void MeanShiftAlgo::filter()开始初始化像素");
	    rootLog.info("void MeanShiftAlgo::filter()初始初始化像素");
		initializePixel(pixels,source);
		subLog.info("void MeanShiftAlgo::filter()初始化像素值成功");
	    rootLog.info("void MeanShiftAlgo::filter()初始化像素值成功");
		vector<MeanPoint*> meanpoints;
		int repeat = 0;
		while(repeat < iterationNums)
		{
			clearMeanPointsContainer(meanpoints);
			for(int i = 0;i < numOfCenters; i++)
			{
				int X = (rand() % (upperIndex[0]-startIndex[0]+1))+1;
				int Y = (rand() % (upperIndex[1]-startIndex[1]+1))+1;
				int Z = (rand() % (upperIndex[2]-startIndex[2]+1))+1;
				itk::Image<signed short,3>::IndexType pixelIndex;
				pixelIndex[0] = X;
				pixelIndex[1] = Y;
				pixelIndex[2] = Z;
				meanpoints.push_back(new MeanPoint(X-startIndex[0],Y-startIndex[1],Z-startIndex[2],source->GetPixel(pixelIndex)));				
			}
			clearAllPoints(allPoints);
			for(int i = 0; i < numOfCenters; i++)
			{
				meanShift(meanpoints[i],allPoints,size[0],size[1],size[2],pixels,pixelFlag);
			}
			subLog.info("void MeanShiftAlgo::filter()meanshift完成，当前迭代次数");
	        rootLog.info("void MeanShiftAlgo::filter()meanshift完成，当前迭代次数");
			//now assign remaining pixels to reature space
			try{
				for(int i = 0; i < size[0]; i++)
				{
					for(int j = 0; j < size[1]; j++)
					{
						for(int k = 0; k < size[2]; k++)
						{
						
							if(!pixelFlag[i][j][k])
							{
								MeanPoint* smp = findSimilarMeans(pixels[i][j][k],allPoints);
								PixelPoint* pixelPoint = new PixelPoint(i,j,k,pixels[i][j][k]);
								allPoints[smp].insert(pixelPoint);
							}else
							{
								pixelFlag[i][j][k] = false;
							}
						}
					}
					cout <<"   "<<i<<"一层一层来 "<<endl;
				}
			}catch(exception&e)
			{
				cout << "跑大循环出错了"<<e.what()<<endl;
			}
			subLog.info("void MeanShiftAlgo::filter()meanshift完成找像素点");
	        rootLog.info("void MeanShiftAlgo::filter()meanshift完成找像素点");
			//update with means centers
			updateMeansPoint(pixels,allPoints);

			repeat++;
		}
		subLog.info("void MeanShiftAlgo::filter()meanshift完成");
	    rootLog.info("void MeanShiftAlgo::filter()meanshift完成");
		vector<SetPixelPoint> needToAssign;
		MapPixel clusterPoints;

		//merge result,remove number of pixels is less than 100
		movePixelPointToNeedReassign(allPoints,needToAssign,clusterPoints);
		moveNeedAssignPointToClusterPoint(needToAssign,clusterPoints);
		subLog.info("void MeanShiftAlgo::filter()合并结果");
	    rootLog.info("void MeanShiftAlgo::filter()合并结果");

		//set pixel value
		setPixelItkIamge(source,clusterPoints,startIndex);
		clearAllPoints(clusterPoints);

		
		
	   
	}

	virtual ~MeanShiftAlgo()
	{
		/**
		MapPixel::iterator cur = allPoints.begin();
		while(cur != allPoints.end())
		{
			 SetPixelPoint pixels = cur->second;
			 SetPixelPoint::iterator pixel_cur = pixels.begin();
			 while(pixel_cur != pixels.end())
			 {
				 delete  (*pixel_cur);
				 pixel_cur++;
			 }
			 delete cur->first;
			 cur++;

		}
		*/
	}

private:	
	void setPixelItkIamge(itk::Image<signed short,3>*dest,MapPixel& res,itk::Image<signed short,3>::IndexType startIndex)
	{
		MapPixel::iterator cur = res.begin();
		while(cur != res.end())
		{
			 SetPixelPoint pixels = cur->second;
			 SetPixelPoint::iterator pixel_cur = pixels.begin();
			 while(pixel_cur != pixels.end())
			 {
				 itk::Image<signed short,3>::IndexType pixelIndex;
				 pixelIndex[0] = (*pixel_cur)->X+startIndex[0];
				 pixelIndex[1] = (*pixel_cur)->Y+startIndex[1];
				 pixelIndex[2] = (*pixel_cur)->Z+startIndex[2];
				 dest->SetPixel(pixelIndex,(*pixel_cur)->gray);
				 pixel_cur++;
			 }
			
			 cur++;

		}	

	}



	void moveNeedAssignPointToClusterPoint(vector<SetPixelPoint>&assign,MapPixel&res)
	{
		for(int i = 0; i < assign.size(); i++)
		{
			SetPixelPoint::iterator cur = assign[i].begin();
			while(cur != assign[i].end())
			{
				MeanPoint* meanPoint = findSimilarMeans((*cur)->gray,res);
				if(meanPoint != 0)
				{
					res[meanPoint].insert(*cur);
				}else
				{

				}
			}
		}
	}
	void movePixelPointToNeedReassign(MapPixel&all,vector<SetPixelPoint>&assign,MapPixel&res)
	{
		MapPixel::iterator cur = all.begin();
		while(cur != all.end())
		{
		
			if(cur->second.size() > clusterPointSize)
			{
				res[cur->first] = cur->second;
			}else
			{
				assign.push_back(cur->second);
			}
			cur++;

		}	
		
	}
	void updateMeansPoint(vector<vector<vector<float>>>&gridPixel,MapPixel&points)
	{
		MapPixel::iterator cur = points.begin();
		while(cur != points.end())
		{
			 SetPixelPoint pixels = cur->second;
			 SetPixelPoint::iterator pixel_cur = pixels.begin();
			 while(pixel_cur != pixels.end())
			 {
				 int X = (*pixel_cur)->X;
				 int Y = (*pixel_cur)->Y;
				 int Z = (*pixel_cur)->Z;
				 float gray = cur->first->gray;
				 gridPixel[X][Y][Z] = gray;
				 pixel_cur++;
			 }
			
			 cur++;

		}	
		
	}

	MeanPoint* findSimilarMeans(int gray,MapPixel& points)
	{
		if(points.size() == 0)
		{
		    return 0;		  
		}
		MapPixel::iterator cursor = points.begin();
		double minDis = DBL_MAX;
		MeanPoint* res = 0;
		while(cursor != points.end())
		{
			double deltgray = gray - cursor->first->gray;
			double dis = deltgray * deltgray;
			if(dis < minDis)
			{
				res = cursor->first;
				minDis = dis;
			}
			cursor++;
		}
		return res;
	}

	//这里传进来的meanPoint什么时候释放
	void meanShift(MeanPoint* meanPoint,MapPixel& points,int XLength,int YLength,int ZLength,
		vector<vector<vector<float>>>&pixels,vector<vector<vector<bool>>>&marked)
	{
		float shift = 0.0;

		//space distance and color distance
		float radius2 = radius * radius;
		float dis2 = grayDistance * grayDistance;

		//current shift center
		int xc = meanPoint->X;
		int yc = meanPoint->Y;
		int zc = meanPoint->Z;
		int grayc = meanPoint->gray;

		//save previous center
		int xcold = meanPoint->X;
		int ycold = meanPoint->Y;
		int zcold = meanPoint->Z;

		float grayOld = meanPoint->gray;

		SetPixelPoint results;

		do
		{
			xcold = xc;
			ycold = yc;
			zcold = zc;
			grayc = grayOld;

			float mx = 0;
			float my = 0;
			float mz = 0;
			float mgray = 0;

			int num = 0;
			// calculate the sum based on generated pixel
			for(int rx = - radius; rx <= radius; rx++)
			{
				int x2 = xc + rx;
				if(x2 >= 0 && x2 < XLength)
				{
					for(int ry = -radius; ry <= radius; ry++)
					{
						int y2 = yc + ry;
						if(y2 >= 0 && y2 < YLength)
						{
							for(int rz = -radius; rz <= radius; rz++)
							{
								int z2 = zc + rz;
								if( z2 >= 0 && z2 < ZLength)
								{
									if(rx*rx + ry*ry +rz*rz <= radius2)
									{
										float gray2 = pixels[x2][y2][z2];
										float dgray = grayc - gray2;
										if(dgray*dgray <= dis2)
										{
											PixelPoint* f = new PixelPoint(x2,y2,z2,gray2);
											results.insert(f);
											marked[x2][y2][z2] = true;
											mx += x2;
											my += y2;
											mz += z2;
											mgray += gray2;
											num++;

										}


									}
								}
							}
						}
					}

				}
			}

			float num_ = 1.0 / num;
			grayc = mgray * num_;
			xc = (int)(mx*num + 0.5);
			yc = (int)(my*num + 0.5);
			zc = (int)(mz*num + 0.5);

			int dx = xc - xcold;
			int dy = yc - ycold;
			int dz = zc - zcold;
			int dgray = grayc - grayOld;

			//shift
			shift = dx * dx + dy * dy + dz * dz;

			//update center location and gray value
			meanPoint->gray = grayc;
			meanPoint->X = xc;
			meanPoint->Y = yc;
			meanPoint->Z = zc;

		} while (shift>0.1);

		//start to merge the features,  find the local maximum
		bool flag = false;
		MapPixel::iterator cur = points.begin();
		while(cur != points.end())
		{
			int deltaX = cur->first->X - meanPoint->X;
			int deltaY = cur->first->Y - meanPoint->Y;
			int deltaZ = cur->first->Z - meanPoint->Z;

			float deltaGray = cur->first->gray - meanPoint->gray;
			float twoSpaceDis = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;
			float twoGrayDis = deltaGray * deltaGray;
			if(twoSpaceDis < radius2 && twoGrayDis < dis2)
			{
				mergeSetPixelPointIntoFirst(cur->second,results);
				flag = true;
				break;
			}

		}
		if(!flag)
		{
			points[meanPoint] = results;
		}



	}

	//inAllPoints属于所有点里某个meanPoint的像素点，meanShiftResult针对某个meanPoint跑出来的结果
	void mergeSetPixelPointIntoFirst(SetPixelPoint&inAllPoints,SetPixelPoint&meanShiftResult)
	{
		set<PixelPoint*>::iterator cur = meanShiftResult.begin();
		while(cur != meanShiftResult.end())
		{
			if(inAllPoints.find(*cur) != inAllPoints.end())
			{
				delete (*cur);
			}else
			{
				inAllPoints.insert(*cur);
			}
		}
	}
	void clearAllPoints(MapPixel&target)
	{
		MapPixel::iterator cur = target.begin();
		while(cur != target.end())
		{
			 SetPixelPoint pixels = cur->second;
			 SetPixelPoint::iterator pixel_cur = pixels.begin();
			 while(pixel_cur != pixels.end())
			 {
				 delete  (*pixel_cur);
				 pixel_cur++;
			 }
			 delete cur->first;
			 cur++;

		}	
		target.clear();
	}
	void clearMeanPointsContainer(vector<MeanPoint*> container)
	{
		for(int i = 0; i < container.size(); i++)
		{
			if(container[i] != 0)
			{
				delete container[i];
			}
		}
		container.clear();
	}
	void initializePixel(vector<vector<vector<float>>>&pixel,itk::Image<signed short,3>*source)
	{
			itk::Image<signed short,3>::IndexType startIndex = source->GetLargestPossibleRegion().GetIndex();
		    itk::Image<signed short,3>::IndexType upperIndex = source->GetLargestPossibleRegion().GetUpperIndex();			
			for(int i = startIndex[0]; i <=upperIndex[0]; i++)
			{
				
				for(int j = startIndex[1]; j <= upperIndex[1]; j++)
				{
					for(int k = startIndex[2]; k <= upperIndex[2]; k++)
					{
						//cout<< i <<" " <<j << " "<< k<<endl;
						itk::Image<signed short,3>::IndexType pixelIndex;
						pixelIndex[0] = i;
						pixelIndex[1] = j;
						pixelIndex[2] = k;
						//subLog.info("void initializePixel(vector<vector<vector<float>>>&pixel,itk::Image<signed short,3>*source)获得某点像素值"+toString(source->GetPixel(pixelIndex)));
						pixel[i-startIndex[0]][j-startIndex[1]][k-startIndex[2]] = source->GetPixel(pixelIndex);
					}
				}
			}
			
		
	}
	




};



#endif