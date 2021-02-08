#include "Simulation.h"
#include <stdio.h>
#include <omp.h>


Simulation::Simulation(Domain& domain)
	: m_domain(domain)
{
}

Simulation::~Simulation()
{
}

void Simulation::Iterate()
{
	Macroscopic();
	Collision();
	PeriodicBC();
	Streaming();
}

void Simulation::Collision()
{
#pragma omp parallel for
	for (int i = 0; i < m_domain.m_domainSize[0] * m_domain.m_domainSize[1]; i++)
	{
		for (int dir = 0; dir < m_domain.GetVelSet()->Get_nDirections(); dir++)
		{
			m_domain.m_lattice[i]->m_newDistributions[dir] = m_domain.m_lattice[i]->m_distributions[dir] - 1 / m_relaxation * (m_domain.m_lattice[i]->m_distributions[dir] - m_domain.m_lattice[i]->Equilibrium(m_domain.GetVelSet(), dir, m_force)) 
				+ m_domain.m_lattice[i]->Force(m_domain.GetVelSet(), dir, m_force);		// applying BGK approximation
		}
	}
}

void Simulation::Streaming()
{
#pragma omp parallel for
	for (int i = 0; i < m_domain.m_domainSize[0] * m_domain.m_domainSize[1]; i++)
	{
		for (int dir = 1; dir < m_domain.GetVelSet()->Get_nDirections(); dir++)
		{
			if (m_domain.m_lattice[i]->m_neighbours[dir] != nullptr)
				m_domain.m_lattice[i]->Stream(dir);
			else if (m_domain.m_lattice[i]->m_neighbours[dir] == nullptr && m_domain.m_lattice[i]->x_position != 0 && m_domain.m_lattice[i]->x_position != m_domain.m_domainSize[0] - 1)	//walls
			{
				int opp_dir = m_domain.GetVelSet()->OppositeDirection(dir);
				m_domain.m_lattice[i]->BounceBack(dir, opp_dir);
			}
		}
	}
}

void Simulation::PeriodicBC()
{
#pragma omp parallel for
	for (int coord_y = 0; coord_y < m_domain.m_domainSize[1]; coord_y++)
	{
		int idx_in = coord_y * m_domain.m_domainSize[0];
		int idx_out = (coord_y + 1) * m_domain.m_domainSize[0] - 1;

		//2nd order
		/*for (int dir = 1; dir < m_domain.GetVelSet()->Get_nDirections(); dir++)
		{
			double w = m_domain.m_velSet->GetWeight(dir);
			double cx = m_domain.m_velSet->GetDirection(dir)[0];
			double cy = m_domain.m_velSet->GetDirection(dir)[1];
			double ux_out = m_domain.m_lattice[idx_out - 1]->m_velocity[0];
			double uy_out = m_domain.m_lattice[idx_out - 1]->m_velocity[1];
			double ux_in = m_domain.m_lattice[idx_in + 1]->m_velocity[0];
			double uy_in = m_domain.m_lattice[idx_in + 1]->m_velocity[1];
			m_domain.m_lattice[idx_in]->m_distributions[dir] = w * (1 + 3 * (cx*ux_out + cy * uy_out) + m_domain.m_lattice[idx_out]->m_distributions[dir] - m_domain.m_lattice[idx_out]->Equilibrium(m_domain.m_velSet, dir, m_force));
			m_domain.m_lattice[idx_out]->m_distributions[dir] = w * (1 + 3 * (cx*ux_in + cy * uy_in) + m_domain.m_lattice[idx_in]->m_distributions[dir] - m_domain.m_lattice[idx_in]->Equilibrium(m_domain.m_velSet, dir, m_force));
		}*/

		//1st order
		m_domain.m_lattice[idx_in]->m_distributions[1] = m_domain.m_lattice[idx_out]->m_distributions[1];
		m_domain.m_lattice[idx_in]->m_distributions[5] = m_domain.m_lattice[idx_out]->m_distributions[5];
		m_domain.m_lattice[idx_in]->m_distributions[8] = m_domain.m_lattice[idx_out]->m_distributions[8];
		m_domain.m_lattice[idx_out]->m_distributions[3] = m_domain.m_lattice[idx_in]->m_distributions[3];
		m_domain.m_lattice[idx_out]->m_distributions[6] = m_domain.m_lattice[idx_in]->m_distributions[6];
		m_domain.m_lattice[idx_out]->m_distributions[7] = m_domain.m_lattice[idx_in]->m_distributions[7];
	}
}

void Simulation::Macroscopic()
{
#pragma omp parallel for
	for (int i = 0; i < m_domain.m_domainSize[0] * m_domain.m_domainSize[1]; i++)
	{
		m_domain.m_lattice[i]->Density();
		m_domain.m_lattice[i]->Velocity(m_domain.m_velSet, m_force);
	}
}

void Simulation::SetRelaxation(double tau)
{
	m_relaxation = tau;
}

void Simulation::SetForce(double force)
{
	m_force = force;
	m_domain.m_force = force;
}

void Simulation::WriteOutput(int iter, double u_max)
{
	OutputFile file((iter));
	file.WriteCSV(m_domain.m_lattice, m_domain.m_velSet, u_max, m_domain.m_domainSize, m_force);
}