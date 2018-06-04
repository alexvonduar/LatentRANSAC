#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif



#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>

#include "config/ConfigParams.h"
#include "estimators/FundmatrixEstimator.h"
#include "estimators/HomogEstimator.h"
#include "estimators/PnPEstimator.h"
#include "estimators/Rigid3dEstimator.h"

// helper functions
bool readCorrsFromFile(std::string& inputFilePath, std::vector<double>& pointData, unsigned int& numPts, bool is3D)
{
	// read data from from file
	std::ifstream infile(inputFilePath.c_str());
	if (!infile.is_open())
	{
		std::cerr << "Error opening input points file: " << inputFilePath << std::endl;
		return false;
	}
	infile >> numPts;
	pointData.resize(6 * numPts);
	for (unsigned int i = 0; i < numPts; ++i)
	{
		if (is3D)
			infile >> pointData[6 * i] >> pointData[6 * i + 1] >> pointData[6 * i + 2] >> pointData[6 * i + 3] >> pointData[6 * i + 4] >> pointData[6 * i + 5];
		else // 2D
		{
			infile >> pointData[6 * i] >> pointData[6 * i + 1] >> pointData[6 * i + 3] >> pointData[6 * i + 4];
			pointData[6 * i + 2] = 1.0;
			pointData[6 * i + 5] = 1.0;
		}
	}
	infile.close();
	return true;
}

bool readPROSACDataFromFile(std::string& sortedPointsFile, unsigned int numPts, std::vector<unsigned int>& prosacData)
{
	std::ifstream infile(sortedPointsFile.c_str());
	if (!infile.is_open())
	{
		std::cerr << "Error opening sorted indices file: " << sortedPointsFile << std::endl;
		return false;
	}
	for (unsigned int i = 0; i < numPts; ++i)
	{
		infile >> prosacData[i];
	}
	infile.close();
	return true;
}


// ---------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	// check command line args
	if (argc < 3)
	{
		std::cerr << "Usage: RunSingleTest <estimation problem> <config file>" << std::endl;
		std::cerr << "\t<estimation problem>: 0 (fundamental matrix), 1 (homography), 2 (pnp), 3 (rigid3d)" << std::endl;
		std::cerr << "\t<config file>: full path to configuration file" << std::endl;
		return(EXIT_FAILURE);
	}
	int estimation_problem = atoi(argv[1]);
	std::string cfg_file_path = argv[2];

	// seed random number generator
	srand((unsigned int)time(NULL));

	if (estimation_problem == 0)
		std::cout << "------- Solving FUNDAMENTAL --------" << std::endl;
	else if (estimation_problem == 1)
		std::cout << "------- Solving HOMOGRAPHY --------" << std::endl;
	else if (estimation_problem == 2)
		std::cout << "------- Solving PNP --------" << std::endl;
	else if (estimation_problem == 3)
		std::cout << "------- Solving rigid-3D --------" << std::endl;
	else
	{
		std::cerr << "Error - not an implemented estimation problem" << std::endl;
		return(EXIT_FAILURE);
	}

	// initialize the appropriate robust estimation problem
	if (estimation_problem == 0)  // FUNDAMENTAL //
	{
		// ------------------------------------------------------------------------
		// initialize the fundamental matrix estimation problem
		ConfigParamsFund cfg;
		if ( !cfg.initParamsFromConfigFile(cfg_file_path) )
		{
			std::cerr << "Error during initialization" << std::endl;
			return(EXIT_FAILURE);
		}
		FundMatrixEstimator* fund = new FundMatrixEstimator;
		fund->initParamsUSAC(cfg);

		// read in input data points
		bool is3D = false; // 2D pixels
		std::vector<double> point_data;
		if (!readCorrsFromFile(cfg.fund.inputFilePath, point_data, cfg.common.numDataPoints, is3D))
		{
			return(EXIT_FAILURE);
		}

		// read in prosac data if required
		std::vector<unsigned int> prosac_data;
		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
			prosac_data.resize(cfg.common.numDataPoints);
			if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
			{
				return(EXIT_FAILURE);
			}
			cfg.prosac.sortedPointIndices = &prosac_data[0];
		} else {
			cfg.prosac.sortedPointIndices = NULL;
		}

		// set up the fundamental matrix estimation problem

		fund->initDataUSAC(cfg);
		fund->initProblem(cfg, &point_data[0]);
		if (!fund->solve())		
			return(EXIT_FAILURE);

		// write out results
		unsigned int pos = (cfg.fund.inputFilePath).find_last_of("/\\");
		std::string working_dir = (cfg.fund.inputFilePath).substr(0, pos + 1);
		std::ofstream outmodel((working_dir + "F.txt").c_str());
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 3; ++j)
			{
				outmodel << fund->final_model_params_[3*i+j] << " ";
			}
		}
		outmodel.close();
		std::ofstream outinliers((working_dir + "inliers.txt").c_str());
		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
		{
			outinliers << fund->usac_results_.inlier_flags_[i] << std::endl;
		}
		outinliers.close();

		// clean up
		point_data.clear();
		prosac_data.clear();
		fund->cleanupProblem();
		delete fund;

	} else if (estimation_problem == 1) { // HOMOGRAPHY //
		// ------------------------------------------------------------------------
		// initialize the homography estimation problem
		ConfigParamsHomog cfg;
		if ( !cfg.initParamsFromConfigFile(cfg_file_path) )
		{
			std::cerr << "Error during initialization" << std::endl;
			return(EXIT_FAILURE);
		}

		HomogEstimator* homog = new HomogEstimator;
		homog->initParamsUSAC(cfg);

		// read in input data points
		bool is3D = false; // 2D pixels
		std::vector<double> point_data;
		if (!readCorrsFromFile(cfg.homog.inputFilePath, point_data, cfg.common.numDataPoints, is3D))
		{
			return(EXIT_FAILURE);
		}

		// read in prosac data if required
		std::vector<unsigned int> prosac_data;
		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
			prosac_data.resize(cfg.common.numDataPoints);
			if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
			{
				return(EXIT_FAILURE);
			}
			cfg.prosac.sortedPointIndices = &prosac_data[0];
		} else {
			cfg.prosac.sortedPointIndices = NULL;
		}

		// set up the homography estimation problem
		homog->initDataUSAC(cfg);
		homog->initProblem(cfg, &point_data[0]);
		if (!homog->solve())
			return(EXIT_FAILURE);

		// write out results
		unsigned int pos = (cfg.homog.inputFilePath).find_last_of("/\\");
		std::string working_dir = (cfg.homog.inputFilePath).substr(0, pos + 1);
		std::ofstream outmodel((working_dir + "H.txt").c_str());
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 3; ++j)
			{
				outmodel << homog->final_model_params_[3*i+j] << " ";
			}
		}
		outmodel.close();
		std::ofstream outinliers((working_dir + "inliers.txt").c_str());
		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
		{
			outinliers << homog->usac_results_.inlier_flags_[i] << std::endl;
		}
		outinliers.close();

		// clean up
		point_data.clear();
		prosac_data.clear();
		homog->cleanupProblem();
		delete homog;
	}	else if (estimation_problem == 2) { // PnP //
		// ------------------------------------------------------------------------
		// initialize the PnP estimation problem
		ConfigParamsPnP cfg;
		if (!cfg.initParamsFromConfigFile(cfg_file_path))
		{
			std::cerr << "Error during initialization" << std::endl;
			return(EXIT_FAILURE);
		}

		PnPEstimator* pnp = new PnPEstimator;
		pnp->initParamsUSAC(cfg);

		// read in input data points
		bool is3D = true; // 3d points and 2d given in homogeneous coordinates
		std::vector<double> point_data;
		if (!readCorrsFromFile(cfg.pnp.inputFilePath, point_data, cfg.common.numDataPoints, is3D))
		{
			return(EXIT_FAILURE);
		}

		// read in prosac data if required
		std::vector<unsigned int> prosac_data;
		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
			prosac_data.resize(cfg.common.numDataPoints);
			if (!readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data))
			{
				return(EXIT_FAILURE);
			}
			cfg.prosac.sortedPointIndices = &prosac_data[0];
		}
		else {
			cfg.prosac.sortedPointIndices = NULL;
		}

		// set up the pnp estimation problem
		pnp->initDataUSAC(cfg);
		pnp->initProblem(cfg, &point_data[0]);

		if (!pnp->solve())
			return(EXIT_FAILURE);

		// write out results
		unsigned int pos = (cfg.pnp.inputFilePath).find_last_of("/\\");
		std::string working_dir = (cfg.pnp.inputFilePath).substr(0, pos + 1);
		std::ofstream outmodel((working_dir + "RT.txt").c_str());
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 4; ++j)
			{
				outmodel << pnp->final_model_params_[i + 3*j] << " "; // save in transposed form
			}
		}
		outmodel.close();
		std::ofstream outinliers((working_dir + "inliers.txt").c_str());
		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
		{
			outinliers << pnp->usac_results_.inlier_flags_[i] << std::endl;
		}
		outinliers.close();

		// clean up
		point_data.clear();
		prosac_data.clear();
		pnp->cleanupProblem();
		delete pnp;

	}
	else if (estimation_problem == 3) { // Rigid 3D //
										// ------------------------------------------------------------------------
										// initialize the Rigid 3D estimation problem
		ConfigParamsRigid3d cfg;
		if (!cfg.initParamsFromConfigFile(cfg_file_path))
		{
			std::cerr << "Error during initialization" << std::endl;
			return(EXIT_FAILURE);
		}

		Rigid3dEstimator* rigid3d = new Rigid3dEstimator;
		rigid3d->initParamsUSAC(cfg);

		// read in input data points
		bool is3D = true; // 3d points
		std::vector<double> point_data;
		if (!readCorrsFromFile(cfg.rigid3d.inputFilePath, point_data, cfg.common.numDataPoints, is3D))
			return(EXIT_FAILURE);

		// read in prosac data if required
		std::vector<unsigned int> prosac_data;
		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
			prosac_data.resize(cfg.common.numDataPoints);
			if (!readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data))
			{
				return(EXIT_FAILURE);
			}
			cfg.prosac.sortedPointIndices = &prosac_data[0];
		}
		else {
			cfg.prosac.sortedPointIndices = NULL;
		}

		// set up the rigid3d estimation problem
		rigid3d->initDataUSAC(cfg);
		rigid3d->initProblem(cfg, &point_data[0]);

		if (!rigid3d->solve())
			return(EXIT_FAILURE);

		// write out results
		unsigned int pos = (cfg.rigid3d.inputFilePath).find_last_of("/\\");
		std::string working_dir = (cfg.rigid3d.inputFilePath).substr(0, pos + 1);
		std::ofstream outmodel((working_dir + "RT.txt").c_str());
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 4; ++j)
			{
				outmodel << rigid3d->final_model_params_[i + 3 * j] << " "; // save in transposed form
			}
		}
		outmodel.close();
		std::ofstream outinliers((working_dir + "inliers.txt").c_str());
		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
		{
			outinliers << rigid3d->usac_results_.inlier_flags_[i] << std::endl;
		}
		outinliers.close();

		// clean up
		point_data.clear();
		prosac_data.clear();
		rigid3d->cleanupProblem();
		delete rigid3d;
	}
	else {
		std::cout << "Estimation problem currently not implemented" << std::endl;
	}

	return(EXIT_SUCCESS);
}

