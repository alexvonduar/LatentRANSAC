#ifndef HOMOGESTIMATOR_H
#define HOMOGESTIMATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <assert.h> 
#include "config/ConfigParamsHomog.h"
#include "utils/MathFunctions.h"
#include "utils/FundmatrixFunctions.h"
#include "utils/HomographyFunctions.h"
#include "USAC.h"

inline bool isClockWise(double x1, double y1, double x2, double y2, double x3, double y3)
{
	return (y3 - y2)*(x2 - x1) < (y2 - y1)*(x3 - x2);
}

class HomogEstimator: public USAC<HomogEstimator>
{
	public:
		inline bool initProblem(const ConfigParamsHomog& cfg, double* pointData);
		// ------------------------------------------------------------------------
		// storage for the final result
		std::vector<double> final_model_params_;

	public:
		HomogEstimator() 
		{
			input_points_ = NULL;
			data_matrix_  = NULL;
			models_.clear();
			models_denorm_.clear();
			// problem dependent params from base class
			hash_dim_ = 8; // hashing dimension
			model_size_ = 9; // number of elements in model matrix

		};
		~HomogEstimator() 
		{
			if (input_points_) { delete[] input_points_; input_points_ = NULL; }
			if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
			for (unsigned int i = 0; i < models_.size(); ++i)
			{
				if (models_[i]) { delete[] models_[i]; }
			}
			models_.clear();
			for (unsigned int i = 0; i < models_denorm_.size(); ++i)
			{
				if (models_denorm_[i]) { delete[] models_denorm_[i]; }
			}
			models_denorm_.clear();
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

	private:
		double*      input_points_denorm_;					// stores pointer to original input points

		// ------------------------------------------------------------------------
		// temporary storage
		double* input_points_;							// stores normalized data points
		double* data_matrix_;							// linearized input data
		double  m_T1_[9], m_T2_[9], m_T2inv_[9];			// normalization matrices
		std::vector<double*> models_;				    // stores vector of models
		std::vector<double*> models_denorm_;			// stores vector of (denormalized) models
		double minminPosition_; // used as base location for hashing
		double minPosition2D_x_, minPosition2D_y_, maxPosition2D_x_, maxPosition2D_y_; // ranges for target quad (also used in hashing)
};


// ============================================================================================
// initProblem: initializes problem specific data and parameters
// this function is called once per run on new data
// ============================================================================================
bool HomogEstimator::initProblem(const ConfigParamsHomog& cfg, double* pointData)
{

	// target ranges - Perhaps should be moved somewhere else
	double exceedingFact = 0.8;
	minPosition2D_x_ = -exceedingFact*latent_image_width_;
	minPosition2D_y_ = -exceedingFact*latent_image_height_;
	maxPosition2D_x_ = (1 + exceedingFact)*latent_image_width_;
	maxPosition2D_y_ = (1 + exceedingFact)*latent_image_height_;
	minminPosition_ = std::min(minPosition2D_x_, minPosition2D_y_);
	double maxRange = std::max(maxPosition2D_x_ - minPosition2D_x_, maxPosition2D_y_ - minPosition2D_y_);
	int bitsPerCoord = ceil(log2(maxRange / latent_bin_size_));
	assert(bitsPerCoord <= 8); // so we can use uint8

	// copy pointer to input data
	input_points_denorm_ = pointData;
	input_points_       = new double[6*cfg.common.numDataPoints];
	if (input_points_denorm_ == NULL)
	{
		std::cerr << "Input point data not properly initialized" << std::endl;
		return false;
	}
	if (input_points_ == NULL)
	{
		std::cerr << "Could not allocate storage for normalized data points" << std::endl;
		return false;
	}

	// normalize input data
	// following this, input_points_ has the normalized points and input_points_denorm_ has 
	// the original input points
	FTools::normalizePoints(input_points_denorm_, input_points_, cfg.common.numDataPoints, m_T1_, m_T2_);
	for (unsigned int i = 0; i < 9; ++i)
	{
		m_T2inv_[i] = m_T2_[i];
	}
	MathTools::minv(m_T2inv_, 3);

	// allocate storage for models
	final_model_params_.clear(); final_model_params_.resize(9);
	models_.clear(); models_.resize(usac_max_solns_per_sample_);
	models_denorm_.clear(); models_denorm_.resize(usac_max_solns_per_sample_);
	for (unsigned int i = 0; i < usac_max_solns_per_sample_; ++i)
	{
		models_[i] = new double[9];
		models_denorm_[i] = new double[9];
	}

	// precompute the data matrix
	data_matrix_ = new double[18*usac_num_data_points_];	// 2 equations per correspondence
	HTools::computeDataMatrix(data_matrix_, usac_num_data_points_, input_points_);

	return true;
}


// ============================================================================================
// cleanupProblem: release any temporary problem specific data storage 
// this function is called at the end of each run on new data
// ============================================================================================
void HomogEstimator::cleanupProblem()
{
	if (input_points_) { delete[] input_points_; input_points_ = NULL; }
	if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
	for (unsigned int i = 0; i < models_.size(); ++i)
	{
		if (models_[i]) { delete[] models_[i]; }
	}
	models_.clear();
	for (unsigned int i = 0; i < models_denorm_.size(); ++i)
	{
		if (models_denorm_[i]) { delete[] models_denorm_[i]; }
	}
	models_denorm_.clear();
}


// ============================================================================================
// generateMinimalSampleModels: generates minimum sample model from the data points whose  
// indices are currently stored in min_sample_. 
// in this case, only one model per minimum sample
// ============================================================================================
unsigned int HomogEstimator::generateMinimalSampleModels()
{
   double A[8*9];
   double At[9*8];

	// form the matrix of equations for this minimal sample
	double *src_ptr;
	double *dst_ptr = A;
	for (unsigned int i = 0; i < usac_min_sample_size_; ++i)
	{
		for (unsigned int j = 0; j < 2; ++j)
		{
			src_ptr = data_matrix_ + 2*min_sample_[i] + j;
			for (unsigned int k = 0; k < 9; ++k)
			{
				*dst_ptr = *src_ptr; 
				++dst_ptr;
				src_ptr += 2*usac_num_data_points_;
			}
		}
	}

	MathTools::mattr(At, A, 8, 9);

	double D[9], U[9*9], V[8*8], *p;
	MathTools::svduv(D, At, U, 9, V, 8);
	p = U + 8;

	double T2_H[9];
	for (unsigned int i = 0; i < 9; ++i)
	{
		*(models_[0]+i) = *p;
		p += 9;
	}
	MathTools::mmul(T2_H, m_T2inv_, models_[0], 3);
	MathTools::mmul(models_denorm_[0], T2_H, m_T1_, 3);  

	return 1;
}


// ============================================================================================
// hashMinimalSampleModels: inserts current sample model(s) into hash  
// 
// returns number of collisions
// ============================================================================================
unsigned int HomogEstimator::hashMinimalSampleModels(unsigned int modelIndex)
{
	double* model = models_denorm_[modelIndex];
	// TODO - remember to create some hash tables upon initialization

	int n_collisions = 0;

	// generate corners of "image1"
	double im_size[2];
	im_size[0] = latent_image_width_;	im_size[1] = latent_image_height_;

	double im_corners[8];
	//im_corners[0] = 1;				im_corners[1] = 1;
	//im_corners[2] = 1;				im_corners[3] = im_size[1];
	//im_corners[4] = im_size[0];		im_corners[5] = im_size[1];
	//im_corners[6] = im_size[0];		im_corners[7] = 1;

	im_corners[0] = 0.2*im_size[0];		im_corners[1] = 0.2*im_size[1];
	im_corners[2] = 0.2*im_size[0];		im_corners[3] = 0.8*im_size[1];
	im_corners[4] = 0.8*im_size[0];		im_corners[5] = 0.8*im_size[1];
	im_corners[6] = 0.8*im_size[0];		im_corners[7] = 0.2*im_size[1];

	// generate quad by transforming corners of "image1"
	double res_quad[8];
	double im_cor[3], vmul_res[3];
	im_cor[2] = 1;
	for (int c = 0; c < 4; c++)
	{
		im_cor[0] = im_corners[2 * c]; im_cor[1] = im_corners[2 * c + 1];
		MathTools::vmul(vmul_res, model, im_cor, 3);
		for (int i = 0; i < 2; i++)
		{
			res_quad[2 * c + i] = vmul_res[i] / vmul_res[2]; // im_size[i] / 2;
			//cout << res_quad[2 * c + i] << " ";
		}
		//cout << endl;
	}

	// validate that quad is not skewed or 'rocket (v) shaped'
	//xL = xL(:, [1:4 1 2]);
	int cycInds[6] = {0,1,2,3,0,1};
	bool areCW;
	for (int c = 0; c < 4; c++)
	{
		areCW = isClockWise(
			res_quad[2 * cycInds[c]], res_quad[2 * cycInds[c] + 1],
			res_quad[2 * cycInds[c + 1]], res_quad[2 * cycInds[c + 1] + 1],
			res_quad[2 * cycInds[c + 2]], res_quad[2 * cycInds[c + 2] + 1]);

		if (!areCW)
			return 0;

	}

	// check range
	bool inRange = true;
	for (int jj = 0; jj < 4; jj++)
	{
		inRange = inRange && res_quad[2 * jj] > minPosition2D_x_ && res_quad[2 * jj] < maxPosition2D_x_;
		inRange = inRange && res_quad[2 * jj + 1] > minPosition2D_y_ && res_quad[2 * jj + 1] < maxPosition2D_y_;
	}

	// such a target location is not valid and we cannot hash it anyway (as we won't have valid values in [0,255])
	if (!inRange)
		return 0;

	for (unsigned int i = 0; i < 2 * usac_min_sample_size_; i++)
		res_quad[i] -= minminPosition_;


	// store hypotesis in history
	//int cur_hyp_index_ = latent_hyp_history_.size();
	//LatentHypothesis h;
	//h.allocate(usac_min_sample_size_, hash_dim_, 9);//(num samples, quad size, model size)

	cur_hyp_index_++;

	unsigned int *ptr1 = &(h_sample_indices_[usac_min_sample_size_ * cur_hyp_index_]);
	for (unsigned int i = 0; i < usac_min_sample_size_; i++)
	{
		//h.sample_indices_[i] = min_sample_[i];
		ptr1[i] = min_sample_[i];
	}

	double *ptr2 = &(h_model_[9 * cur_hyp_index_]);
	for (int i = 0; i < 9; i++)
	{
		//h.model_[i] = model[i];
		ptr2[i] = model[i];
	}

	ptr2 = &(h_hash_seed_[hash_dim_ * cur_hyp_index_]);
	for (int i = 0; i < hash_dim_; i++)
	{
		//h.hash_seed_[i] = res_quad[i];
		ptr2[i] = res_quad[i];
	}

	for (int i = 0; i < 8; i++)
		res_quad[i] = res_quad[i] / latent_bin_size_;
	//latent_hyp_history_.push_back(h);


	// quantize resulting quad & insert into rand-grid(s)
	for (unsigned int g = 0; g < latent_num_grids_; g++)
	{
		// quantize resulting quad
		unsigned long long hash_64 = 0;
		unsigned char* ptr8 = (unsigned char *)(&hash_64);
		int ptr = hash_dim_ * g;

		for (int i = 0; i < 8; i++)
		{
			*ptr8 = (unsigned char)(res_quad[i] + latent_grid_shifts_[ptr++]);
			ptr8++;
		}

		// calc hash value
		//unsigned int hash_val = hash64toX(hash_64,latent_num_bits_);
		unsigned int hash_val = hash64byMod(hash_64, latent_num_bits_);

		// push into hash table & check for collision
		int jjj;
		int tableWidth = 50;
		for (jjj = 0; jjj < tableWidth; jjj++)
		{

			int collision_index = latent_rand_grids_[g][hash_val + jjj];
			if (collision_index != -1)
			{ // -> this index was already "hit"  before - validate it is a real collision

				// validate if this is really a collision
				if (validateCollision(collision_index, cur_hyp_index_))
				{
					//cur_hyp_collisions.push_back(collision_index);
					n_collisions++;
					// update parent/child indices
					//latent_hyp_history_[cur_hyp_index].father_ = collision_index;
					//latent_hyp_history_[collision_index].children_indices_.push_back(cur_hyp_index);
					break; // don't go through complete table "width"
				}
			}
			else
				break;
		}		
		
		if (jjj < tableWidth)
		{
			// store model in table
			latent_rand_grids_[g][hash_val + std::min(jjj, tableWidth - 1)] = cur_hyp_index_; // mark this location with this hypothesis
		}
		else
		{ // clear and start all over
			for (int kkk = 0; kkk < tableWidth; kkk++)
			latent_rand_grids_[g][hash_val + kkk] = -1;// mark this location with this hypothesis
			// store model in table
			latent_rand_grids_[g][hash_val] = cur_hyp_index_;// mark this location with this hypothesis
		}						

		if (n_collisions)
			break; // don't go through all tables

	}// end for over all grids

	return n_collisions;
}

// ============================================================================================
// validateCollision: check if a collision was caused by two similar hypothesis (and not false)
// i.e. - check if original (seed) matches are different *enough*, and if quads are similar *enough*
// ============================================================================================
inline bool HomogEstimator::validateCollision(int index1, int index2)
{

	unsigned int* h_sample_indices_1 = &(h_sample_indices_[usac_min_sample_size_ * index1]);
	unsigned int* h_sample_indices_2 = &(h_sample_indices_[usac_min_sample_size_ * index2]);


#if 0
	// Bidirectional validity test - check each of the colliding models on the others' seed points
	// Simon: changed here the thresh to rely on the inlier rate threshold 
	// old: double bidirectionalThresh = latent_bidirectional_factor_ * latent_collision_tolerance_;
	double bidirectionalThresh = latent_bidirectional_factor_ * usac_inlier_threshold_square_;
	bool bidirectionalValidity = true;
	//bidirectionalValidity = bidirectionalValidity && checkSampsonVsThreshold(h1.model_, h2.sample_indices_, bidirectionalThresh);
	//bidirectionalValidity = bidirectionalValidity && checkSampsonVsThreshold(h2.model_, h1.sample_indices_, bidirectionalThresh);
	double *h_model_1 = &(h_model_[9 * index1]);
	double *h_model_2 = &(h_model_[9 * index2]);
	bidirectionalValidity = bidirectionalValidity && checkSampsonVsThreshold(h_model_1, h_sample_indices_2, bidirectionalThresh);
	bidirectionalValidity = bidirectionalValidity && checkSampsonVsThreshold(h_model_2, h_sample_indices_1, bidirectionalThresh);
	
	if (!bidirectionalValidity)
		return false;
#endif


#if 1
	double *hash_seed_1 = &(h_hash_seed_[hash_dim_ * index1]);
	double *hash_seed_2 = &(h_hash_seed_[hash_dim_ * index2]);

	for (int i = 0; i < hash_dim_; i++)
	{
		double diff_dist = abs(hash_seed_1[i] - hash_seed_2[i]);
		if (diff_dist > latent_collision_tolerance_)
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

	int maxNumCommon = 2;
	if (v_intersection.size() > maxNumCommon)
		return false;
#endif


	return true;
}

// ============================================================================================
// generateRefinedModel: compute model using non-minimal set of samples
// default operation is to use a weight of 1 for every data point
// ============================================================================================
bool HomogEstimator::generateRefinedModel(std::vector<unsigned int>& sample,
										  unsigned int numPoints,
										  bool weighted,
										  double* weights)
{
	// form the matrix of equations for this non-minimal sample
	double *A = new double[numPoints*2*9];	
	double *src_ptr;
	double *dst_ptr = A;
	for (unsigned int i = 0; i < numPoints; ++i)
	{
		for (unsigned int j = 0; j < 2; ++j)
		{
			src_ptr = data_matrix_ + 2*sample[i] + j;
			for (unsigned int k = 0; k < 9; ++k)
			{
				if (!weighted)
				{
					*dst_ptr = *src_ptr;
				}
				else
				{
					*dst_ptr = (*src_ptr)*weights[i];
				}
				++dst_ptr;
				src_ptr += 2*usac_num_data_points_;
			}
		}
	}

	// decompose
	double V[9*9], D[9], *p;
	MathTools::svdu1v(D, A, 2*numPoints, V, 9);

	unsigned int j = 0;
	for (unsigned int i = 1; i < 9; ++i)
	{
		if (D[i] < D[j]) 
			j = i;
	}
	p = V + j;

	for (unsigned int i = 0; i < 9; ++i)
	{
		*(models_[0]+i) = *p;
		p += 9;
	}
	double T2_H[9];
	MathTools::mmul(T2_H, m_T2inv_, models_[0], 3);
	MathTools::mmul(models_denorm_[0], T2_H, m_T1_, 3); 

	delete[] A;

	return true;
}


// ============================================================================================
// validateSample: check if minimal sample is valid
// ============================================================================================
bool HomogEstimator::validateSample()
{
	// check oriented constraints
   double p[3], q[3];
   double *a, *b, *c, *d;

   a = input_points_ + 6*min_sample_[0];
   b = input_points_ + 6*min_sample_[1];
   c = input_points_ + 6*min_sample_[2];
   d = input_points_ + 6*min_sample_[3];

   HTools::crossprod(p, a, b, 1);
   HTools::crossprod(q, a+3, b+3, 1);

   if ((p[0]*c[0]+p[1]*c[1]+p[2]*c[2])*(q[0]*c[3]+q[1]*c[4]+q[2]*c[5])<0)
      return false;
   if ((p[0]*d[0]+p[1]*d[1]+p[2]*d[2])*(q[0]*d[3]+q[1]*d[4]+q[2]*d[5])<0)
      return false;

   HTools::crossprod(p, c, d, 1);
   HTools::crossprod(q, c+3, d+3, 1);

   if ((p[0]*a[0]+p[1]*a[1]+p[2]*a[2])*(q[0]*a[3]+q[1]*a[4]+q[2]*a[5])<0)
      return false;
   if ((p[0]*b[0]+p[1]*b[1]+p[2]*b[2])*(q[0]*b[3]+q[1]*b[4]+q[2]*b[5])<0)
      return false;

   return true;	
}


// ============================================================================================
// validateModel: check if model computed from minimal sample is valid
// ============================================================================================
bool HomogEstimator::validateModel(const unsigned int modelIndex)
{
	return true;
}


inline bool HomogEstimator::checkSampsonVsThreshold(double* model, unsigned int* samplesInds, double threshold)
{
	double inv_model[9];
	double h_x[3], h_inv_xp[3];
	double* pt;
	std::vector<double>::iterator current_err_array = err_ptr_[0];

	for (unsigned int i = 0; i < 9; ++i)
	{
		inv_model[i] = model[i];
	}
	MathTools::minv(inv_model, 3);

	// compute symmetric transfer error
	for (unsigned int i = 0; i < usac_min_sample_size_; i++)
	{
		pt = input_points_denorm_ + 6 * samplesInds[i];
		MathTools::vmul(h_x, model, pt, 3);

		double err1 = 0.0;
		for (unsigned int j = 0; j < 2; ++j)		
			err1 += (h_x[j] / h_x[2] - pt[3 + j]) * (h_x[j] / h_x[2] - pt[3 + j]);
		
		if (err1 > threshold)
			return false;

		MathTools::vmul(h_inv_xp, inv_model, pt + 3, 3);
		double err2 = 0.0;
		for (unsigned int j = 0; j < 2; ++j)
			err2 += (h_inv_xp[j] / h_inv_xp[2] - pt[j]) * (h_inv_xp[j] / h_inv_xp[2] - pt[j]);

		if (err2 > threshold)
			return false;
	}
	return true;

}

// ============================================================================================
// evaluateModel: test model against all/subset of the data points
// ============================================================================================
bool HomogEstimator::evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested)
{
	double* model = models_denorm_[modelIndex];
	double inv_model[9];
	double h_x[3], h_inv_xp[3], temp_err;
	double* pt;
	std::vector<double>::iterator current_err_array = err_ptr_[0];
	bool good_flag = true;
	double lambdaj, lambdaj_1 = 1.0;
	*numInliers = 0;
	*numPointsTested = 0;
	unsigned int pt_index;

	for (unsigned int i = 0; i < 9; ++i)
	{
		inv_model[i] = model[i];
	}
	MathTools::minv(inv_model, 3);
	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		// get index of point to be verified
		if (eval_pool_index_ > usac_num_data_points_-1)
		{
			eval_pool_index_ = 0;
		}
		pt_index = evaluation_pool_[eval_pool_index_];
		++eval_pool_index_;

		// compute symmetric transfer error
		pt = input_points_denorm_ + 6*pt_index;
		MathTools::vmul(h_x, model, pt, 3);
		MathTools::vmul(h_inv_xp, inv_model, pt+3, 3);

		double err1 = 0.0, err2 = 0.0;
		for (unsigned int j = 0; j < 2; ++j)
		{
			err1 += (h_x[j]/h_x[2] - pt[3+j]) * (h_x[j]/h_x[2] - pt[3+j]);
			err2 += (h_inv_xp[j]/h_inv_xp[2] - pt[j]) * (h_inv_xp[j]/h_inv_xp[2] - pt[j]);
		}
		temp_err = 0.5f*(err1 + err2); // NOTE: we changed from sum to average
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
void HomogEstimator::testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel)
{
	*degenerateModel = false;
	*upgradeModel = false;
}

// ============================================================================================
// upgradeDegenerateModel: try to upgrade degenerate model to non-degenerate by sampling from
// the set of outliers to the degenerate model
// ============================================================================================
unsigned int HomogEstimator::upgradeDegenerateModel()
{
	return 0;
}


// ============================================================================================
// findWeights: given model and points, compute weights to be used in local optimization
// ============================================================================================
void HomogEstimator::findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
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
void HomogEstimator::storeModel(const unsigned int modelIndex, unsigned int numInliers)
{
	// save the current model as the best solution so far
	for (unsigned int i = 0; i < 9; ++i)
	{
		final_model_params_[i] = *(models_denorm_[modelIndex]+i);
	}
}

#endif

