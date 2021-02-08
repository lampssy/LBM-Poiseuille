#include "Domain.h"


Domain::Domain(std::shared_ptr<VelocitySet> velSet, const std::vector<int> domainSize)
	: m_velSet(velSet), m_domainSize(domainSize)
{
}


Domain::~Domain()
{
}


void Domain::Initialize()
{
	for (auto node_ptr : m_lattice)
	{
		node_ptr->m_distributions.reserve(m_velSet->Get_nDirections());
		for (int dir = 0; dir < m_velSet->Get_nDirections(); dir++)
		{
			node_ptr->Initialize(m_velSet, dir, m_force);
		}
	}
}

void Domain::ConnectNodes()
{
	for (auto node_ptr : m_lattice)
	{
		node_ptr->m_neighbours[0] = nullptr;		// no neighbour in the 0 direction
		for (int dir = 1; dir < m_velSet->Get_nDirections(); dir++)
		{
			int x_neighbour = node_ptr->x_position + m_velSet->GetDirection(dir)[0];
			int y_neighbour = node_ptr->y_position + m_velSet->GetDirection(dir)[1];
			if (IsInDomain(x_neighbour, y_neighbour))
			{
				int idx_neighbour = y_neighbour * m_domainSize[0] + x_neighbour;
				node_ptr->m_neighbours[dir] = m_lattice[idx_neighbour];
			}
			else
				node_ptr->m_neighbours[dir] = nullptr;
		}
	}
}

void Domain::CreateNodes()
{
	int max_nNodes = m_domainSize[1] * m_domainSize[0];
	m_lattice.reserve(max_nNodes);
	for (int coord_y = 0; coord_y < m_domainSize[1]; coord_y++)
	{
		for (int coord_x = 0; coord_x < m_domainSize[0]; coord_x++)
		{
			m_lattice.push_back(std::make_shared<Node>(coord_x, coord_y));		// creating pointers to internal and wall nodes
		}
	}
}

bool Domain::IsInDomain(int coord_x, int coord_y)
{
	if (coord_x >= 0 && coord_x < m_domainSize[0] && coord_y >= 0 && coord_y < m_domainSize[1])
		return 1;
	else
		return 0;
}

const std::shared_ptr<VelocitySet> Domain::GetVelSet()
{
	return m_velSet;
}