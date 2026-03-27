#include "Body.h"

Body::Body(double mass, Vector starting_pos, Vector starting_velocity)
{
	this->mass = mass;
	position = starting_pos;
	velocity = starting_velocity;
}
