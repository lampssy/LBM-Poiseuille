#pragma once
#include "VelocitySet.h"
#include "Node.h"
#include <vector>
#include <iostream>

class Domain
{
public:
	int m_nDimensions = m_domainSize.size();
private:
	std::shared_ptr<VelocitySet> m_velSet;
	std::vector<std::shared_ptr<Node>> m_lattice;
	std::vector<int> m_domainSize;
	double m_force;
public:
	Domain(std::shared_ptr<VelocitySet> velSet, const std::vector<int> domainSize);
	~Domain();

	void CreateNodes();
	bool IsInDomain(int coord_x, int coord_y);
	const std::shared_ptr<VelocitySet> GetVelSet();
	void ConnectNodes();
	void Initialize();

	friend class Simulation;
	friend class Node;
};

