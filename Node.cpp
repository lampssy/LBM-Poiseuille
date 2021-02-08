#include "Node.h"

Node::Node(int coord_x, int coord_y)
	:	x_position(coord_x), y_position(coord_y)
{
}


Node::~Node()
{
}

void Node::Stream(int dir)
{
	m_neighbours[dir]->m_distributions[dir] = m_newDistributions[dir];
}

double Node::Equilibrium(std::shared_ptr<VelocitySet> velSet, int dir, double force)
{
	double du = m_velocity[0] * velSet->GetDirection(dir)[0] + m_velocity[1] * velSet->GetDirection(dir)[1];

	// incompressible linear equilibrium
	return velSet->GetWeight(dir) * (m_density + 3 * du);
}

double Node::Force(std::shared_ptr<VelocitySet> velSet, int dir, double force)
{
	double w = velSet->GetWeight(dir);
	double cx = velSet->GetDirection(dir)[0];
	double cy = velSet->GetDirection(dir)[1];
	return  w * 3 * cx * force; //1st order
	//double ux = m_velocity[0];
	//double uy = m_velocity[1];
	//return w * ((3 * (cx - ux) + 9 * (cx*ux + cy*uy) * force) + 3 * (cy - uy)); //2nd order
}

void Node::Initialize(std::shared_ptr<VelocitySet> velSet, int dir, double force)
{
	double df = force * velSet->GetDirection(dir)[0];
	m_distributions.push_back(velSet->GetWeight(dir) - 0.5 * velSet->GetWeight(dir) * 3 * df);
	m_newDistributions.push_back(velSet->GetWeight(dir) - 0.5 * velSet->GetWeight(dir) * 3 * df);
}

void Node::Density()
{
	double density = 0;
	for (int dir = 0; dir < 9; dir++)
	{
		density += this->m_distributions[dir];	// sum of distributions for all discrete velocities (in all directions)
	}
	m_density = density;
}

void Node::Velocity(std::shared_ptr<VelocitySet> velSet, double force)
{
	std::vector<double> velocity;
	velocity.reserve(2);	// 2D case
	velocity.push_back(0);
	velocity.push_back(0);

	for (int dir = 0; dir < 9; dir++)
	{
		for (int dim = 0; dim < 2; dim++)
		{
			velocity[dim] += velSet->GetDirection(dir)[dim] * this->m_distributions[dir];
		}
	}

	//velocity[0] += force / 2; //2nd order

	m_velocity = velocity;
}

void Node::BounceBack(int dir, int opp_dir)
{
	m_distributions[opp_dir] = m_newDistributions[dir];
}
