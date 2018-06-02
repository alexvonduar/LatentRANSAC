#ifndef PNPESTIMATOR_H
#define PNPESTIMATOR_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;
#include <iterator>
//#include <algorithm>
//#include <math.h>       /* sqrt */
#include <assert.h>      
#include "config/ConfigParamsPnP.h"
#include "utils/MathFunctions.h"
#include "utils/PnpFunctions.h"
//#include "utils/FundmatrixFunctions.h"
#include "USAC.h"
#include <omp.h>

//include generic headers for opengv stuff
#include <opengv/types.hpp>
#include <opengv/absolute_pose/modules/main.hpp>
#include <opengv/absolute_pose/modules/Epnp.hpp>
#include <opengv/math/roots.hpp>
//#include <opengv/math/arun.hpp>
#include <Eigen/src/Geometry/AngleAxis.h>

#define PI 3.141592653589793238463

class PnPEstimator: public USAC<PnPEstimator>
{
	public:
		inline bool		 initProblem(const ConfigParamsPnP& cfg, double* pointData);
		// ------------------------------------------------------------------------
		// storage for the final result
		std::vector<double> final_model_params_;

	public:
		PnPEstimator() 
		{
			models_denorm_.clear();
			// problem dependent params from base class
			hash_dim_ = 6; // hashing dimension
			model_size_ = 12; // number of elements in model matrix

		};
		~PnPEstimator() 
		{	
			//if (!models_denorm_.empty())
			//{
			//	for (unsigned int i = 0; i < models_denorm_.size(); ++i)
			//	{
			//		if (models_denorm_[i]) { delete[] models_denorm_[i]; }
			//	}
			//	models_denorm_.clear();			
			//}

		};

	public:
		// ------------------------------------------------------------------------
		// problem specific functions
		inline void		 cleanupProblem();
		inline unsigned int generateMinimalSampleModels();
		inline bool			validateCollision(int index1, int index2);
		inline unsigned int hashMinimalSampleModels(unsigned int modelIndex);
		inline bool		 generateRefinedModel(std::vector<unsigned int>& sample, const unsigned int numPoints, 
										  bool weighted = false, double* weights = NULL);
		inline bool		 validateSample();
		inline bool		 validateModel(unsigned int modelIndex);
		inline bool		 evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested);
		inline void		 testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel);
		inline unsigned int upgradeDegenerateModel();
		inline void		 findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
									 unsigned int numInliers, double* weights);
		inline void		 storeModel(unsigned int modelIndex, unsigned int numInliers);
		inline bool      checkSampsonVsThreshold(double* model, unsigned int* samplesInds, double threshold);

		inline bool      getPoint3D(unsigned int index, opengv::point_t& point);
		inline bool      getPoint2D(unsigned int index, opengv::point_t& point);

		inline void PnPEstimator::NormalizePoint(double* point);

private:
		inline double getResidualError(double* correspondence, double* R_prime, double* T);
		double*      input_points_denorm_;					// stores pointer to original input points
#if 1
		
		opengv::rotations_t rotations;
		opengv::translations_t translations;
		opengv::transformations_t solutions;
#endif

		// ------------------------------------------------------------------------
		// temporary storage
		std::vector<double*> models_denorm_;			// stores vector of (denormalized) models

		double minPositionAA, minPositionT, maxPositionAA, maxPositionT; // ranges for target latent_vector (also used in hashing)
		double bin_size_AA, bin_size_T; // bin sizes for hashing
		double collision_tolerance_R, collision_tolerance_T; // collision tolerance is different for R and T
};

// ============================================================================================
// NormalizePoint: Normalize a 3D point (a,b,c) by dividing it by sqrt(a^2+b^2+c^2)
// ============================================================================================
inline void PnPEstimator::NormalizePoint(double* point)
{
	double denom = sqrt(point[0] * point[0] + point[1] * point[1] + point[2] * point[2]);
	for (int i = 0; i < 3; i++)
		point[i] = point[i] / denom;
}

inline bool PnPEstimator::getPoint3D(unsigned int index, opengv::point_t& point)
{

	double* basePoint = input_points_denorm_ + 6 * index;

	point[0] = basePoint[0];
	point[1] = basePoint[1];
	point[2] = basePoint[2];
	return true;
}


inline bool PnPEstimator::getPoint2D(unsigned int index, opengv::point_t& point)
{

	double* basePoint = input_points_denorm_ + 6 * index;

	point[0] = basePoint[3];
	point[1] = basePoint[4];
	point[2] = basePoint[5];
	return true;
}

// ============================================================================================
// initProblem: initializes problem specific data and parameters
// this function is called once per run on new data
// ============================================================================================
bool PnPEstimator::initProblem(const ConfigParamsPnP& cfg, double* pointData)
{

	// target ranges - Perhaps should be moved somewhere else
//	double exceedingFact = 0.8;
	minPositionAA = -PI; 
	minPositionT = -latent_maximum_parallax_; 
	maxPositionAA = PI;
	maxPositionT = latent_maximum_parallax_; 

	collision_tolerance_R = cfg.latent.collision_tolerance; 
	collision_tolerance_T = collision_tolerance_R * 2.1; // tan(collision_tolerance_R) * latent_maximum_parallax_;
	cout << "!!!!! collision tolerance ration (R/T) is set to 2.1" << endl;

	bin_size_AA = collision_tolerance_R * cfg.latent.bin_size_factor;
	bin_size_T = collision_tolerance_T * cfg.latent.bin_size_factor;

	cout << "+++++ collision_tolerance_R: " << collision_tolerance_R << "  collision_tolerance_T: " << collision_tolerance_T << endl;
	double range;
	int bitsPerCoord;
	range = maxPositionAA - minPositionAA;
	bitsPerCoord = ceil(log2(range / bin_size_AA));
	assert(bitsPerCoord <= 8); // so we can use uint8
	range = maxPositionT - minPositionT;
	bitsPerCoord = ceil(log2(range / bin_size_T));
	assert(bitsPerCoord <= 8); // so we can use uint8
	usac_inlier_threshold_square_ = sqrt(usac_inlier_threshold_square_); // This is overridden here for PnP - do not take square (since errors aren't square)
	assert(lo_num_iterative_steps_ == 0); // this should be set through the config. We don't handle for pnp weighted lo-ransac

	// copy pointer to input data
	input_points_denorm_ = pointData;
	if (input_points_denorm_ == NULL)
	{
		std::cerr << "Input point data not properly initialized" << std::endl;
		return false;
	}

	// allocate storage for models
	final_model_params_.clear(); final_model_params_.resize(model_size_);
	models_denorm_.clear(); models_denorm_.resize(usac_max_solns_per_sample_);
	for (unsigned int i = 0; i < usac_max_solns_per_sample_; ++i)
	{
		models_denorm_[i] = new double[model_size_]; // R and t
	}

	//// precompute the data matrix

	return true;
}


// ============================================================================================
// cleanupProblem: release any temporary problem specific data storage 
// this function is called at the end of each run on new data
// ============================================================================================
void PnPEstimator::cleanupProblem()
{
	
	for (unsigned int i = 0; i < models_denorm_.size(); ++i)
	{
		if (models_denorm_[i]) { delete[] models_denorm_[i]; }
	}
	models_denorm_.clear();
}


// ============================================================================================
// generateMinimalSampleModels: generates minimum sample model from the data points whose  
// indices are currently stored in min_sample_. 
// computed models are stored in models_denorm_
// retuns the number of fitted models (which in this case might be up to 4)
// ============================================================================================
unsigned int PnPEstimator::generateMinimalSampleModels()
{

	// three 3D points
	opengv::bearingVectors_t f;
	opengv::bearingVector_t bv;
	for (int i = 0; i < 3; i++)
	{
		getPoint3D(min_sample_[i], bv);
		f.push_back(bv);
	}

	// three 2D points
	opengv::points_t p;
	opengv::point_t point;
	for (int i = 0; i < 3; i++)
	{
		getPoint2D(min_sample_[i], point);
		p.push_back(point);
	}

	// solving using Kneip's p3p algorithm
	opengv::rotations_t all_rotations;
	opengv::translations_t all_translations;
	opengv::transformations_t all_solutions;

	solutions.clear();
	rotations.clear();
	translations.clear();
	opengv::absolute_pose::modules::p3p_kneip_main(p, f, all_solutions);
	PTools::p3p_kneip_RT(p, f, all_rotations, all_translations, all_solutions);

	// extract solutions
	int p3p_num_solutions = 4;
	unsigned int nsols = 0; 
	double* result_data;
	double tmp[12];
	double tmp2[12];
	for (int i = 0; i < all_solutions.size(); i++)
	{
		result_data = all_solutions[i].data(); // the .data() gives the transposed matrix, which is actually what we need for later on
		// check for bad (NaN) solution
		if (isnan(result_data[0]) || isnan(result_data[11])) // should be enough to check these two entries of the matrix
			continue;

		memcpy((void*)models_denorm_[nsols], result_data, model_size_ * sizeof(double));
		
		rotations.push_back(all_rotations[i]);
		translations.push_back(all_translations[i]);
		solutions.push_back(all_solutions[i]);

		nsols++;
	}

	// print solutions
	bool printThis = 0;
	if (printThis)
	{
		cout << "Kneip's P3P algorithm produced " << nsols << " actual solutions:" << std::endl;
		for (unsigned int i = 0; i < nsols; i++)
		{
			cout << "--- solution " << i << " ---" << endl;
			cout << solutions[i] << std::endl << std::endl;
			cout << rotations[i] << std::endl << std::endl;
			cout << translations[i] << std::endl << std::endl;
		}
	}
	return nsols;

}


// ============================================================================================
// hashMinimalSampleModels: inserts current sample model(s) into hash  
// 
// returns number of collisions
// ============================================================================================
unsigned int PnPEstimator::hashMinimalSampleModels(unsigned int modelIndex)
{
	int n_collisions = 0;

		// the model to hash
		double* model = models_denorm_[modelIndex];

		opengv::rotation_t R = rotations[modelIndex];
		opengv::translation_t T = translations[modelIndex];

		Eigen::AngleAxis<double> eigen_aa = Eigen::AngleAxisd(R);

		// create an axis-angle representation (make axis have length of angle)
		Eigen::Vector3d AA = eigen_aa.axis() * eigen_aa.angle();

		// centering
		for (int i = 0; i < 3; i++)
		{
			T[i] = T[i] - minPositionT;
			AA[i] = AA[i] - minPositionAA;
		}

	    cur_hyp_index_++;
		// store the sample indices
		unsigned int *ptr1 = &(h_sample_indices_[usac_min_sample_size_ * cur_hyp_index_]);
		for (int i = 0; i < usac_min_sample_size_; i++)
			ptr1[i] = min_sample_[i];

		// store the model (matrix)
		double *ptr2 = &(h_model_[model_size_ * cur_hyp_index_]);
		for (int i = 0; i < model_size_; i++)
			ptr2[i] = model[i];

		// store the hash seed (the vector before quantization, but after 'centering')
		ptr2 = &(h_hash_seed_[hash_dim_ * cur_hyp_index_]);
		int ptr = 0;
		for (int i = 0; i < 3; i++)
			ptr2[ptr++] = AA[i];
		for (int i = 0; i < 3; i++)
			ptr2[ptr++] = T[i]; // including centering


		// binning
		for (int i = 0; i < 3; i++)
		{
			T[i] = T[i] / bin_size_T;
			AA[i] = AA[i] / bin_size_AA;
		}

		// a single latent_vector
		double latent_vector[6];
		for (int i = 0; i < 3; i++)
		{
			latent_vector[i] = AA[i];
			latent_vector[3+i] = T[i];
		}

		// quantize resulting seed & insert into rand-grid(s)
		unsigned long long hash_64 = 0;
		double val_to_quant;
		unsigned int hash_val;
		unsigned int g;
		int i, collision_indices[128]; // 128 is at least the number of grids
		for (g = 0; g < latent_num_grids_; g++)
		{

			// TODO - check if aa vector has magnitude >pi & needs to be "flipped" modulo
			// TODO - check if aa vector has magnitude >pi & needs to be "flipped" modulo
			// TODO - check if aa vector has magnitude >pi & needs to be "flipped" modulo

			// quantize resulting latent_vector and convert to a single large unsigned integer
			hash_64 = 0;
			ptr = hash_dim_ * g;
			unsigned char* ptr8 = (unsigned char *)(&hash_64);
			for (i = 0; i < 6; i++)
			{
				*ptr8 = (unsigned char)(latent_vector[i] + latent_grid_shifts_[ptr++]);
				ptr8++;
				//hash_64 = hash_64 << 8;
				//hash_64 += (static_cast <unsigned long long>(latent_vector[i] + latent_grid_shifts_[ptr++])); // keep only 8 LSB of "val_to_quant" into "hash_64"
			}

			// calc hash value
			hash_val = hash64byMod(hash_64, latent_num_bits_);
			//hash_val = (unsigned int)(hash_64 % 524287);

			// check for a collision ("hit")
			int collision_index = latent_rand_grids_[g][hash_val];
			//collision_indices[g] = latent_rand_grids_[g][hash_val];

			// validate collisions
			if (collision_index != -1)
			{ // -> this index was already "hit"  before - validate it is a real collision
				// validate if this is really a collision
				if (validateCollision(collision_index, cur_hyp_index_))
				{
					n_collisions += 1;
					//cout << "$$ (" << usac_results_.hyp_count_ << ") $$ found collision (BREAKING!!!) at table " << g << " between current index " << cur_hyp_index_ << " and previous index " << collision_index << endl;
					break;
				}

			}			

			// mark the location in the hash table
			latent_rand_grids_[g][hash_val] = cur_hyp_index_;// mark this location with this hypothesis	

		}// end for over all grids

	return n_collisions;

}

// ============================================================================================
// validateCollision: check if a collision was caused by two similar hypothesis (and not false)
// i.e. - check if original (seed) matches are different *enough*, and if quads are similar *enough*
// ============================================================================================
inline bool PnPEstimator::validateCollision(int index1, int index2)
{

	unsigned int* h_sample_indices_1 = &(h_sample_indices_[usac_min_sample_size_ * index1]);
	unsigned int* h_sample_indices_2 = &(h_sample_indices_[usac_min_sample_size_ * index2]);


#if 1
	// Bidirectional validity test - check each of the colliding models on the others' seed points
	// Simon: changed here the thresh to rely on the inlier rate threshold 
	// old: double bidirectionalThresh = latent_bidirectional_factor_ * latent_collision_tolerance_;
	double bidirectionalThresh = latent_bidirectional_factor_ * usac_inlier_threshold_square_;
	bool bidirectionalValidity = true;
	double *h_model_1 = &(h_model_[model_size_ * index1]);
	double *h_model_2 = &(h_model_[model_size_ * index2]);
	bidirectionalValidity = bidirectionalValidity && checkSampsonVsThreshold(h_model_1, h_sample_indices_2, bidirectionalThresh);
	bidirectionalValidity = bidirectionalValidity && checkSampsonVsThreshold(h_model_2, h_sample_indices_1, bidirectionalThresh);
	
	if (!bidirectionalValidity)
		return false;
#endif


#if 1
	double *hash_seed_1 = &(h_hash_seed_[hash_dim_ * index1]);
	double *hash_seed_2 = &(h_hash_seed_[hash_dim_ * index2]);

	// rotation 
	for (int i = 0; i < hash_dim_ / 2; i++)
	{
		double diff_dist = abs(hash_seed_1[i] - hash_seed_2[i]);
		if (diff_dist > collision_tolerance_R)
			return false; //current corner is too far
	}
	// translation
	for (int i = (hash_dim_ / 2) + 1; i < hash_dim_; i++)
	{
		double diff_dist = abs(hash_seed_1[i] - hash_seed_2[i]);
		if (diff_dist > collision_tolerance_T)
			return false; //current corner is too far
	}
#endif


#if 1

	std::vector<unsigned int> indices1;
	indices1.assign(h_sample_indices_1, h_sample_indices_1 + usac_min_sample_size_);
	std::sort(indices1.begin(), indices1.end());

	std::vector<unsigned int> indices2;
	indices2.assign(h_sample_indices_2, h_sample_indices_2 + usac_min_sample_size_);
	std::sort(indices2.begin(), indices2.end());

	std::vector<unsigned int> v_intersection;
	std::set_intersection(indices1.begin(), indices1.end(), indices2.begin(), indices2.end(), std::back_inserter(v_intersection));

	//if (v_intersection.size() != v_intersection2.size())
	//	cout << "PROBLEM" << endl;

	int maxNumCommon = 1;
	if (v_intersection.size() > maxNumCommon)
		return false;
#endif

	return true;
}

// ============================================================================================
// generateRefinedModel: compute model using non-minimal set of samples
// default operation is to use a weight of 1 for every data point
// ============================================================================================
bool PnPEstimator::generateRefinedModel(std::vector<unsigned int>& sample,
										  unsigned int numPoints,
										  bool weighted,
										  double* weights)
{	

	// all 3D points
	opengv::bearingVectors_t f;
	opengv::bearingVector_t bv;
	for (int i = 0; i < numPoints; i++)
	{
		getPoint3D(sample[i], bv);
		f.push_back(bv);
	}


	// all 2D points
	opengv::points_t p;
	opengv::point_t point;
	for (int i = 0; i < numPoints; i++)
	{
		getPoint2D(sample[i], point);
		p.push_back(point);
	}

	// solving using epnp
	opengv::rotation_t rotation;
	opengv::translation_t translation;
	opengv::transformation_t solution;

	PTools::epnp_RT(p, f, rotation, translation, solution);

	// extract solution
	double* result_data;
	double tmp[12];
	double tmp2[12];

	result_data = solution.data(); // the .data() gives the transposed matrix, which is actually what we need for later on
	// check for bad (NaN) solution
	assert(!(isnan(result_data[0]) || isnan(result_data[11]))); // should be enough to check these two entries of the matrix
	//		double* targetAddress = ((double*)models_denorm_[nsols]);
	memcpy((void*)models_denorm_[0], result_data, model_size_ * sizeof(double));


	// print solutions
	bool printThis = 0;
	if (printThis)
	{
		cout << endl << "The epnp algorithm produced this solution:" << std::endl;
		cout << solution << std::endl << std::endl;
		cout << rotation << std::endl << std::endl;
		cout << translation << std::endl << std::endl;
	}
	

	return true;
}


// ============================================================================================
// validateSample: check if minimal sample is valid
// ============================================================================================
bool PnPEstimator::validateSample()
{
	//// check oriented constraints
 //  double p[3], q[3];
 //  double *a, *b, *c, *d;

 //  a = input_points_ + 6*min_sample_[0];
 //  b = input_points_ + 6*min_sample_[1];
 //  c = input_points_ + 6*min_sample_[2];
 //  d = input_points_ + 6*min_sample_[3];

 //  HTools::crossprod(p, a, b, 1);
 //  HTools::crossprod(q, a+3, b+3, 1);

 //  if ((p[0]*c[0]+p[1]*c[1]+p[2]*c[2])*(q[0]*c[3]+q[1]*c[4]+q[2]*c[5])<0)
 //     return false;
 //  if ((p[0]*d[0]+p[1]*d[1]+p[2]*d[2])*(q[0]*d[3]+q[1]*d[4]+q[2]*d[5])<0)
 //     return false;

 //  HTools::crossprod(p, c, d, 1);
 //  HTools::crossprod(q, c+3, d+3, 1);

 //  if ((p[0]*a[0]+p[1]*a[1]+p[2]*a[2])*(q[0]*a[3]+q[1]*a[4]+q[2]*a[5])<0)
 //     return false;
 //  if ((p[0]*b[0]+p[1]*b[1]+p[2]*b[2])*(q[0]*b[3]+q[1]*b[4]+q[2]*b[5])<0)
 //     return false;

   return true;	
}


// ============================================================================================
// validateModel: check if model computed from minimal sample is valid
// ============================================================================================
bool PnPEstimator::validateModel(const unsigned int modelIndex)
{
	return true;
}

double PnPEstimator::getResidualError(double* correspondence, double* R_prime, double* T)
{
	double err, vec1[3], vec2[3], res;
	double* pt3d = correspondence;
	double* pt2d = correspondence + 3;

	MathTools::subtractvecs(vec1, pt3d, T, 3);
	MathTools::vmul(vec2, R_prime, vec1, 3);
	NormalizePoint(vec2);
	MathTools::dotprod(res, vec2, pt2d, 3);
	//assert(res > -1.005 && res < 1.005);
	res = std::min(1.0, res);
	res = std::max(-1.0, res);
	err = acos(res);
	return err;
}

inline bool PnPEstimator::checkSampsonVsThreshold(double* model, unsigned int* samplesInds, double threshold)
{
	double vec[3], temp_err;
	double *pt,  *pt3d, *pt2d;

	double* R_prime = model;
	double* T = model + 9;
	double res, err;

	// compute angular error
	for (int i = 0; i < usac_min_sample_size_; i++)
	{
		pt = input_points_denorm_ + 6 * samplesInds[i];

		// in a function:
		err = getResidualError(pt, R_prime, T);

		// compare to threshold
		if (err > threshold)
			return false;

	}
	return true;

}

// ============================================================================================
// evaluateModel: test model against all/subset of the data points
// ============================================================================================
bool PnPEstimator::evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested)
{

	double* model = models_denorm_[modelIndex];
	double* pt;
	std::vector<double>::iterator current_err_array = err_ptr_[0];
	bool good_flag = true;
	double lambdaj, lambdaj_1 = 1.0;
	*numInliers = 0;
	*numPointsTested = 0;
	unsigned int pt_index;

	double temp_err;

	double* R_prime = model;
	double* T = model + 9;

	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		// get index of point to be verified
		if (eval_pool_index_ > usac_num_data_points_-1)
		{
			eval_pool_index_ = 0;
		}
		pt_index = evaluation_pool_[eval_pool_index_];
		++eval_pool_index_;

		pt = input_points_denorm_ + 6*pt_index;

		// compute residual error
		temp_err = getResidualError(pt, R_prime, T);

		*(current_err_array+pt_index) = temp_err;

		if (temp_err < usac_inlier_threshold_square_)
		{
			++(*numInliers);
		}

		if (usac_verif_method_ == USACConfig::VERIF_SPRT)
		{
			if (temp_err < usac_inlier_threshold_square_)
			{			
				lambdaj = lambdaj_1 * (sprt_delta_/sprt_epsilon_);
			}
			else
			{
				lambdaj = lambdaj_1 * ( (1 - sprt_delta_)/(1 - sprt_epsilon_) );
			}

			if (lambdaj > decision_threshold_sprt_)
			{
				good_flag = false;
				*numPointsTested = i+1;
				return good_flag;
			}
			else
			{
				lambdaj_1 = lambdaj;
			}
		}
	}
	*numPointsTested = usac_num_data_points_;

	return good_flag;
}

// ============================================================================================
// testSolutionDegeneracy: check if model is degenerate
// ============================================================================================
void PnPEstimator::testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel)
{
	*degenerateModel = false;
	*upgradeModel = false;
}

// ============================================================================================
// upgradeDegenerateModel: try to upgrade degenerate model to non-degenerate by sampling from
// the set of outliers to the degenerate model
// ============================================================================================
unsigned int PnPEstimator::upgradeDegenerateModel()
{
	return 0;
}


// ============================================================================================
// findWeights: given model and points, compute weights to be used in local optimization
// ============================================================================================
void PnPEstimator::findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
								 unsigned int numInliers, double* weights)
{
	for (unsigned int i = 0; i < numInliers; ++i)
	{
		weights[i] = 1.0;
	}
}


// ============================================================================================
// storeModel: stores current best model
// this function is called  (by USAC) every time a new best model is found
// ============================================================================================
void PnPEstimator::storeModel(const unsigned int modelIndex, unsigned int numInliers)
{
	// save the current model as the best solution so far
	for (unsigned int i = 0; i < model_size_; ++i)
	{
		final_model_params_[i] = *(models_denorm_[modelIndex]+i);
	}
}

#endif

