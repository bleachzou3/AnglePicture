#ifndef MEAN_POINT_HPP_
#define MEAN_POINT_HPP_
#include "Vector3.h"
class MeanPoint
{
public:
  int X;
  int Y;
  int Z;
  float gray;
  MeanPoint(int _X,int _Y,int _Z,float _gray)
  {
	  this->X = _X;
	  this->Y = _Y;
	  this->Z = _Z;
	  this->gray = _gray;
  }
};

class PixelPoint
{
public:
	int X;
	int Y;
	int Z;
	float gray;
	int flag;
	PixelPoint(int _X, int _Y, int _Z,float _gray)
	{
	  this->X = _X;
	  this->Y = _Y;
	  this->Z = _Z;
	  this->gray = _gray;
	  flag = -1;
	}
};
#endif