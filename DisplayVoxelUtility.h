#ifndef DISPLAY_VOXEL_UTILITY_HPP
#define DISPLAY_VOXLE_UTILITY_HPP
using namespace std;
#include <string>
class DisplayVoxelUtility
{
public:
	DisplayVoxelUtility();
	virtual ~DisplayVoxelUtility();
	static void displaySegmentBloodVesselsFromVti(string fileName);
};

#endif