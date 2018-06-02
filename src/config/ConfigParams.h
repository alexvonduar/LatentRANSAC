#ifndef CONFIGPARAMS_H
#define CONFIGPARAMS_H

#include <string>
#include "ConfigFileReader.h"

namespace USACConfig
{
	enum RandomSamplingMethod	  {SAMP_UNIFORM, SAMP_PROSAC, SAMP_UNIFORM_MM};
	enum VerifMethod			  {VERIF_STANDARD, VERIF_SPRT};
	enum LocalOptimizationMethod  {LO_NONE, LO_LOSAC};

	// common USAC parameters
	struct Common
	{
		// change these parameters according to the problem of choice
		Common(): confThreshold		     (0.95), 
			      minSampleSize		     (7),
				  inlierThreshold		 (0.001),
				  maxHypotheses		     (100000),
  				  maxSolutionsPerSample  (3),
				  numDataPoints          (0),
				  prevalidateSample 	 (false),
				  prevalidateModel	     (false),
				  testDegeneracy	 	 (false),
				  randomSamplingMethod   (SAMP_UNIFORM),
				  verifMethod			 (VERIF_STANDARD),
				  localOptMethod         (LO_NONE)
		{}
		double				    confThreshold;
		unsigned int		    minSampleSize;  
		double				    inlierThreshold;   
		unsigned int		    maxHypotheses;	
		unsigned int		    maxSolutionsPerSample; 
		unsigned int		    numDataPoints;
		bool					prevalidateSample;
		bool					prevalidateModel;
		bool					testDegeneracy;
		RandomSamplingMethod    randomSamplingMethod;
		VerifMethod			    verifMethod;
		LocalOptimizationMethod localOptMethod;
	};

	// LATENT parameters
	struct Latent
	{
		Latent() : is_active(false),
					image_width			(0),
					image_height		(0),
					maximum_parallax	(0),
					bin_size_factor		(2.0),
					num_grids			(8),
					max_num_collisions  (10),
					table_size_factor(2.0),
					bidirectional_factor(2.0),
					collision_tolerance(10.0),
					random_grid_prob(0.1),
					stop_crit_factor(3)
		{}
		bool			is_active;
		unsigned int	image_width;
		unsigned int	image_height;
		double			maximum_parallax;
		double			bin_size_factor;
		unsigned int	num_grids;
		unsigned int    max_num_collisions;
		double			table_size_factor;
		double			bidirectional_factor;
		double			collision_tolerance;
		double			random_grid_prob;
		double			stop_crit_factor;
	};


	// PROSAC parameters
	struct Prosac
	{
		Prosac(): maxSamples 		    (200000),
				  beta                  (0.05),
				  nonRandConf           (0.95),
				  minStopLen	        (20),
				  sortedPointsFile		(""),		// leave blank if not reading from file
				  sortedPointIndices    (NULL)		// this should point to an array of point indices
													// sorted in decreasing order of quality scores
		{}
		unsigned int  maxSamples;
		double		  beta;
		double        nonRandConf;
		unsigned int  minStopLen;
		std::string   sortedPointsFile;
		unsigned int* sortedPointIndices;
	};

	// SPRT parameters
	struct Sprt
	{
		Sprt(): tM      (200.0),
				mS	    (2.38),
				delta   (0.05),
				epsilon (0.2)
		{}
		double tM;
		double mS;
		double delta;
		double epsilon;
	};

	// LOSAC parameters
	struct Losac
	{
		Losac(): innerSampleSize		  (14),
		 		 innerRansacRepetitions   (10),
				 thresholdMultiplier	  (2.0),
				 numStepsIterative	      (4)
		{}
		unsigned int innerSampleSize;
		unsigned int innerRansacRepetitions;
		double		 thresholdMultiplier;
		unsigned int numStepsIterative;
	};

}

// main configuration struct that is passed to USAC
class ConfigParams
{
public:
	// to be overridden to read in model specific data
	virtual bool initParamsFromConfigFile(std::string& configFilePath);

	// function to read in common usac parameters from config file
	bool initUSACParamsFromConfigFile(const ConfigFileReader& cfr);

	USACConfig::Common     common;
	USACConfig::Latent     latent;
	USACConfig::Prosac     prosac;
	USACConfig::Sprt       sprt;
	USACConfig::Losac      losac;
};


#endif