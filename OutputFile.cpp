#include "OutputFile.h"


OutputFile::OutputFile(int outputIdx)
	:	m_outputIdx(outputIdx)
{
}


OutputFile::~OutputFile()
{
}


void OutputFile::WriteCSV(std::vector<std::shared_ptr<Node>> lattice, std::shared_ptr<VelocitySet> velSet, double u_max, std::vector<int> m_domainSize, double force)
{
	std::ofstream velocityCSV;
	std::string file_name = std::to_string(m_outputIdx) + ".csv";
	velocityCSV.open(file_name);
	//velocityCSV << "coordX" << "," << "coordY" << "," << "coordZ" << "," << "U" << "," << "p" << "," << "f0" << "," << "f1" << "," << "f2" << "," << "f3" << "," << "f4" << "," << "f5" << "," << "f6" << "," << "f7" << "," << "f8" << std::endl;
	velocityCSV << "coordX" << "," << "coordY" << "," << "coordZ" << "," << "p" << "," << "U" << "," << "U_an" << "," << "error" << std::endl;
	int dy = m_domainSize[1];
	int y_max = dy - 1;

	for (auto node : lattice)
	{
		double velMag = sqrt(node->m_velocity[0] * node->m_velocity[0] + node->m_velocity[1] * node->m_velocity[1]);
		double velAnalytical = -4. * u_max / (dy*dy) * (node->y_position + 0.5 - 0.0) * (node->y_position - 0.5 - y_max);
		double error = abs(velAnalytical - velMag)/velAnalytical;
		//velocityCSV << node->x_position << "," << node->y_position << "," << 0 << "," << velMag << "," << node->m_density << "," << node->m_distributions[0] << "," << node->m_distributions[1] << "," << node->m_distributions[2] << "," 
			//<< node->m_distributions[3] << "," << node->m_distributions[4] << "," << node->m_distributions[5] << "," << node->m_distributions[6] << "," << node->m_distributions[7] << "," << node->m_distributions[8] << std::endl;
		velocityCSV << node->x_position << "," << node->y_position << "," << 0 << "," << node->m_density << "," << velMag << "," << velAnalytical << "," << error << std::endl;
	}
}
