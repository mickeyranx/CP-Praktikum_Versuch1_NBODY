#include "vector.h"
#include <cmath>
using namespace std;

namespace Customvectors {


	Vector::Vector(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->length = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	}

	// Hallo Welt 

	Vector::~Vector()
	{
		x = 0;
		y = 0;
		z = 0;
		length = 0;
	}


	Vector::Vector()
	{
	}


	Vector Vector::operator+(Vector vec) {
		return Vector(x + vec.x, y + vec.y, z + vec.z);
	}

	Vector Vector::operator-(Vector vec)
	{
		return Vector(x - vec.x, y - vec.y, z -  vec.z);
	}

	double Vector::operator*(Vector vec) {
		return x * vec.x + y * vec.y + z * vec.z;
	}

	Vector Vector::operator*(double s)
	{
		return Vector(x*s, y*s, z*s);
	}

	Vector Vector::operator%(Vector vec) {
		return Vector(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x);
	}


}
