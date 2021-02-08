#pragma once

#include "VelocitySet.h"
#include "Domain.h"
#include "OutputFile.h"

class Simulation
{
private:
	Domain m_domain;
	double m_relaxation;
	double m_force;
private:
	void Collision();
	void Streaming();
	void PeriodicBC();

public:
	Simulation(Domain& domain);
	~Simulation();

	void Iterate();
	void SetRelaxation(double tau);
	void SetForce(double force);
	void WriteOutput(int iter, double u_max);
	void Macroscopic();
};

