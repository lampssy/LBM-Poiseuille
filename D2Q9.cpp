#include "D2Q9.h"

D2Q9::D2Q9()
	: VelocitySet(
		{ 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0 , 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 },
		{ {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1} },
		1.0 / 3.0
	)
{
}


D2Q9::~D2Q9()
{
}

int D2Q9::OppositeDirection(int dir) const
{
	std::vector<int> v_opp = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
	return v_opp[dir];
}