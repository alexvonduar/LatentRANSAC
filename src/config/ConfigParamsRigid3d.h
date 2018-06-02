#ifndef CONFIGPARAMSRIGID3D_H
#define CONFIGPARAMSRIGID3D_H

#include "ConfigParams.h"

namespace USACConfig
{
	// problem specific/data-related parameters: fundamental matrix
	struct Rigid3d
	{
		Rigid3d(): inputFilePath	      ("")			// leave blank if not using config file
		{}

		std::string			inputFilePath;
	};
}

class ConfigParamsRigid3d: public ConfigParams
{
public:
	// simple function to read in parameters from config file
	bool initParamsFromConfigFile(std::string& configFilePath);

	USACConfig::Rigid3d rigid3d;
};

#endif