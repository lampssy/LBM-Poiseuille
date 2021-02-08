#pragma once
#include <iostream>
#include "VelocitySet.h"

class Node
{
public:
	std::vector<double> m_distributions;
	std::vector<double> m_newDistributions;
	std::shared_ptr<Node> m_neighbours[9];
	double m_density;
	std::vector<double> m_velocity;
	int x_position;
	int y_position;
public:
	Node(int coord_x, int coord_y);
	~Node();

	double Equilibrium(std::shared_ptr<VelocitySet> velSet, int dir, double force);
	double Force(std::shared_ptr<VelocitySet> velSet, int dir, double force);
	void Density();
	void Velocity(std::shared_ptr<VelocitySet> velSet, double force);
	void Stream(int dir);
	void BounceBack(int dir, int opp_dir);
	void Initialize(std::shared_ptr<VelocitySet> velSet, int dir, double force);
};

