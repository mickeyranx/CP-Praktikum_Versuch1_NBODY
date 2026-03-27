#pragma once
/*
this is a custom vector class
for making life easier
*/
namespace Customvectors{


	class Vector
	{
	private:
		double length;
		double x;
		double y;
		double z;
	public:
		Vector(double x, double y, double z);
		~Vector();
		Vector();
		double getX() { return x; };
		double getY() { return y; };
		double getZ() { return z; };
		void setX(double x) { this->x = x; };
		void setY(double y) { this->y = y; };
		void setZ(double z) { this->z = z; };
		double getLength() { return length; };
		
		//Vector addition
		Vector operator+(Vector v);

		//Vector subtraction
		Vector operator-(Vector v);

		//TODO: overload for multiplication with scalar
		//Vector multiplication
		double operator*(Vector vec);

		Vector operator*(double s);

		//Vector product
		Vector operator%(Vector vec);







	};

}