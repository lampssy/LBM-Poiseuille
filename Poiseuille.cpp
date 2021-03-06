// LBM-2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "D2Q9.h"
#include "Simulation.h"
#include "OutputFile.h"

int main()
{
	int scale = 5;
	int dx = 20 * scale;
	int dy = 5 * scale;
	double tau = sqrt(3./16.) + 0.5;
	double u_max = 0.1 / scale;
	double nu = (2 * tau - 1) / 6.;
	double Re = dy * u_max / nu;
	double Fx = 8 * nu*u_max / double(dy*dy);

	std::vector<int> domainSize{ dx, dy };
	std::shared_ptr<D2Q9> velSet = std::make_shared<D2Q9>();
	Domain domain(velSet, domainSize);
	domain.CreateNodes();
	domain.ConnectNodes();

	Simulation sim(domain);
	sim.SetRelaxation(tau);
	sim.SetForce(Fx);

	domain.Initialize();

	int max_iter = 0;
	//std::cout << "Enter the number of iterations: " << std::endl;
	//std::cin >> max_iter;
	max_iter = 10000;

	for (int iter = 1; iter <= max_iter; iter++)
	{
		sim.Iterate();
		std::cout << "Iteration: " << iter << std::endl;

		if (iter % 1000 == 0)
			sim.WriteOutput(iter, u_max);
	}

	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
