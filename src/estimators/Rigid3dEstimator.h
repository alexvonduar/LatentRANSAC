#ifndef RIGID3DESTIMATOR_H
#define RIGID3DESTIMATOR_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;
//#include <algorithm>
//#include <math.h>       /* sqrt */
#include <assert.h>
#include "config/ConfigParamsRigid3d.h"
#include "utils/MathFunctions.h"
#include "utils/FundmatrixFunctions.h"
#include "USAC.h"
#include <omp.h>

//include generic headers for opengv stuff
#include <opengv/types.hpp>
#include <opengv/absolute_pose/modules/main.hpp>
#include <opengv/math/roots.hpp>
#include <opengv/math/arun.hpp>
#include <Eigen/src/Geometry/AngleAxis.h>

#define PI 3.141592653589793238463


#define SQR(a) ((a)*(a))


class Rigid3dEstimator: public USAC<Rigid3dEstimator>
{
    public:
        inline bool initProblem(const ConfigParamsRigid3d& cfg, double* pointData);
        // ------------------------------------------------------------------------
        // storage for the final result
        std::vector<double> final_model_params_;

    public:
        Rigid3dEstimator()
        {
            models_denorm_.clear();
            // problem dependent params from base class
            hash_dim_ = 6; // hashing dimension
            model_size_ = 12; // number of elements in model matrix
            usac_max_solns_per_sample_ = 1;
        };
        ~Rigid3dEstimator()
        {
            //if (!models_denorm_.empty())
            //{
            // for (unsigned int i = 0; i < models_denorm_.size(); ++i)
            // {
            //     if (models_denorm_[i]) { delete[] models_denorm_[i]; }
            // }
            // models_denorm_.clear();
            //}

        };

    public:
        // ------------------------------------------------------------------------
        // problem specific functions
        inline void cleanupProblem();
        inline unsigned int generateMinimalSampleModels();
        inline bool validateCollision(int index1, int index2);
        inline unsigned int hashMinimalSampleModels(unsigned int modelIndex);
        inline bool generateRefinedModel(std::vector<unsigned int>& sample, const unsigned int numPoints,
                                          bool weighted = false, double* weights = NULL);
        inline bool validateSample();
        inline bool validateModel(unsigned int modelIndex);
        inline bool evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested);
        inline void testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel);
        inline unsigned int upgradeDegenerateModel();
        inline void findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers,
                                     unsigned int numInliers, double* weights);
        inline void storeModel(unsigned int modelIndex, unsigned int numInliers);
        inline bool checkSampsonVsThreshold(double* model, unsigned int* samplesInds, double threshold);

        inline void getPointPair3D(unsigned int index, opengv::point_t& point1, opengv::point_t& point2);

        inline void NormalizePoint(double* point);

    private:
        inline double DistR3(double* point1, double* point2);
        inline double getResidualError(double* correspondence, double* R_prime, double* T);
        double*      input_points_denorm_; // stores pointer to original input points

        // ------------------------------------------------------------------------
        // temporary storage
        std::vector<double*> models_denorm_; // stores vector of (denormalized) models

        opengv::rotation_t rotation_;
        opengv::translation_t translation_;


        double minPositionR, minPositionT, maxPositionR, maxPositionT; // ranges for target latent_vector (also used in hashing)
        double congruency_threshold_ = 0.9f; // used for sample validation
        double bin_size_R_, bin_size_T_; // bin sizes for hashing
        double collision_tolerance_R_, collision_tolerance_T_; // collision tolerance is different for R and T
};


// ============================================================================================
// DistR3: Euclidean distance between two points in 3D
// ============================================================================================
inline double Rigid3dEstimator::DistR3(double* point1, double* point2) {
    return sqrt(
        SQR(point1[0] - point2[0]) +
        SQR(point1[1] - point2[1]) +
        SQR(point1[2] - point2[2]));
}


// ============================================================================================
// NormalizePoint: Normalize a 3D point (a,b,c) by dividing it by sqrt(a^2+b^2+c^2)
// ============================================================================================
inline void Rigid3dEstimator::NormalizePoint(double* point)
{
    double denom = sqrt(point[0] * point[0] + point[1] * point[1] + point[2] * point[2]);
    for (int i = 0; i < 3; i++)
        point[i] = point[i] / denom;
}

inline void Rigid3dEstimator::getPointPair3D(unsigned int index, opengv::point_t& point1, opengv::point_t& point2)
{

    double* basePoint = input_points_denorm_ + 6 * index;

    point1[0] = basePoint[0];
    point1[1] = basePoint[1];
    point1[2] = basePoint[2];
    point2[0] = basePoint[3];
    point2[1] = basePoint[4];
    point2[2] = basePoint[5];
    return;
}



// ============================================================================================
// initProblem: initializes problem specific data and parameters
// this function is called once per run on new data
// ============================================================================================
bool Rigid3dEstimator::initProblem(const ConfigParamsRigid3d& cfg, double* pointData)
{

    // target ranges - Perhaps should be moved somewhere else
    // double exceedingFact = 0.8;
    minPositionR = -PI; // -exceedingFact*latent_image_width_;
    minPositionT = -latent_maximum_parallax_; // -exceedingFact*latent_image_height_;
    maxPositionR = PI; // (1 + exceedingFact)*latent_image_width_;
    maxPositionT = latent_maximum_parallax_; // (1 + exceedingFact)*latent_image_height_;
    // minminPosition_ = min(minPositionR, minPositionT);

    collision_tolerance_T_ = cfg.latent.collision_tolerance;
    collision_tolerance_R_ = collision_tolerance_T_ * 3.6f;  // empirically validated on synthetic data
    cout << "!!!!! collision tolerance ration (T/R) is set to 3.6" << endl;

    bin_size_R_ = collision_tolerance_R_ * cfg.latent.bin_size_factor;
    bin_size_T_ = collision_tolerance_T_ * cfg.latent.bin_size_factor;

    cout << "+++++ collision_tolerance_R_: " << collision_tolerance_R_ << "  collision_tolerance_T_: " << collision_tolerance_T_ << endl;
    double range;
    int bitsPerCoord;
    range = maxPositionR - minPositionR;
    bitsPerCoord = ceil(log2(range / bin_size_R_));
    assert(bitsPerCoord <= 8); // so we can use uint8
    range = maxPositionT - minPositionT;
    bitsPerCoord = ceil(log2(range / bin_size_T_));
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
void Rigid3dEstimator::cleanupProblem()
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
unsigned int Rigid3dEstimator::generateMinimalSampleModels()
{

    // copy three 3D point pairs into eign
    opengv::points_t pts1, pts2;
    for (int i = 0; i < 3; i++)
    {
        opengv::point_t p1, p2;
        getPointPair3D(min_sample_[i], p1, p2);
        pts1.push_back(p1);
        pts2.push_back(p2);
    }


    // solving using arun's algorithm
    opengv::transformation_t solution;
    solution = opengv::math::arun_complete(pts1, pts2);
    rotation_ = solution.block<3, 3>(0, 0);
    translation_ = solution.col(3);

    memcpy((void*)models_denorm_[0], solution.data(), model_size_ * sizeof(double));

    return 1;
}


// ============================================================================================
// hashMinimalSampleModels: inserts current sample model(s) into hash
//
// returns number of collisions
// ============================================================================================
unsigned int Rigid3dEstimator::hashMinimalSampleModels(unsigned int modelIndex)
{
    int n_collisions = 0;


    // the model to hash
    double* model = models_denorm_[modelIndex];

    opengv::rotation_t R = rotation_;
    opengv::translation_t T = translation_;

    Eigen::AngleAxis<double> eigen_aa = Eigen::AngleAxisd(R);

    // create an axis-angle representation (make axis have length of angle)
    Eigen::Vector3d AA = eigen_aa.axis() * eigen_aa.angle();

    // T centering
    for (int i = 0; i < 3; i++)
        T[i] = T[i] - minPositionT;

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
        ptr2[ptr++] = AA[i] - minPositionR;
    for (int i = 0; i < 3; i++)
        ptr2[ptr++] = T[i]; // including centering


                            // binning
    for (int i = 0; i < 3; i++)
    {
        T[i] = T[i] / bin_size_T_;
        AA[i] = AA[i] / bin_size_R_;
    }

    // a single latent_vector
    double latent_vector[6];
    for (int i = 0; i < 3; i++)
    {
        latent_vector[i] = AA[i];
        latent_vector[3 + i] = T[i];
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
            n_collisions += validateCollision(collision_index, cur_hyp_index_);
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
inline bool Rigid3dEstimator::validateCollision(int index1, int index2)
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
        if (diff_dist > collision_tolerance_R_)
            return false; //current corner is too far
    }
    // translation
    for (int i = (hash_dim_ / 2) + 1; i < hash_dim_; i++)
    {
        double diff_dist = abs(hash_seed_1[i] - hash_seed_2[i]);
        if (diff_dist > collision_tolerance_T_)
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
    //    cout << "PROBLEM" << endl;

    int maxNumCommon = 2;
    if (v_intersection.size() > maxNumCommon)
        return false;
#endif




    //LatentHypothesis h1 = latent_hyp_history_[index1];
    //LatentHypothesis h2 = latent_hyp_history_[index2];


    //// STEP 1: comapre quads, too different quads are false collision
    //for (int i = 0; i < 4; i++)
    //{
    //    double dx = h1.hash_seed_[2 * i] - h2.hash_seed_[2 * i];
    //    double dy = h1.hash_seed_[2 * i + 1] - h2.hash_seed_[2 * i + 1];

    //    if (sqrt(dx*dx + dy*dy) > latent_collision_tolerance_)
    //        return false;//current corner is too far
    //}

    //// STEP 2: check that the two samples of matches don't have too much in common (not just identity, but actual location)
    //// sort needed before using "next_permutation" in the next loop
    //std::sort(h1.sample_indices_.begin(), h1.sample_indices_.end());

    //int numClosePairs = 0;
    //do // run through all permuations
    //{
    //    double absDiffs[8]; // = 2 * usac_min_sample_size_

    //    for (int i = 0; i < usac_min_sample_size_; i++)
    //    {
    //        int i1 = h1.sample_indices_[i];
    //        int i2 = h2.sample_indices_[i];
    //        double* pt1 = input_points_denorm_ + 6 * i1;
    //        double* pt2 = input_points_denorm_ + 6 * i2;

    //        double sqrDiffLeft,sqrDiffRight;

    //        sqrDiffLeft = mySqr(*(pt1 + 0) - *(pt2 + 0));
    //        sqrDiffLeft += mySqr(*(pt1 + 1) - *(pt2 + 1));
    //
    //        sqrDiffRight = mySqr(*(pt1 + 3) - *(pt2 + 3));
    //        sqrDiffRight += mySqr(*(pt1 + 4) - *(pt2 + 4));

    //        numClosePairs += (sqrDiffLeft < usac_inlier_threshold_square_) && (sqrDiffRight < usac_inlier_threshold_square_);
    //    }

    //    if (numClosePairs > 1) // more than 1 out of 4 close points
    //        return false;//at least two matches are close to identical
    //} while (std::next_permutation(h1.sample_indices_.begin(), h1.sample_indices_.end()));

#ifdef DUMP_FOR_DEBUG
    // FOUND A COLLISSION - DUMP SOME DEBUG INFO


    std::ofstream ofs1("c:\\temp\\matrix1TMP.txt");
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 4; ++j)
        {
            ofs1 << h_model_1[4 * i + j] << " ";
        }
    }
    ofs1.close();

    std::ofstream ofs2("c:\\temp\\matrix2TMP.txt");
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 4; ++j)
        {
            ofs2 << h_model_2[4 * i + j] << " ";
        }
    }
    ofs2.close();

    std::ofstream ofs3("c:\\temp\\LRc_inds1TMP.txt");
    for (int i = 0; i < usac_min_sample_size_; i++)
        ofs3 << h_sample_indices_1[i] << std::endl;
    ofs3.close();

    std::ofstream ofs4("c:\\temp\\LRc_inds2TMP.txt");
    for (int i = 0; i < usac_min_sample_size_; i++)
        ofs4 << h_sample_indices_2[i] << std::endl;
    ofs4.close();

    std::ofstream ofs5("c:\\temp\\collisionIndsTMP.txt");
    ofs5 << index1 << std::endl << index2 << std::endl;
    ofs5.close();

    std::ofstream ofs6("c:\\temp\\LRc_quad1TMP.txt");
    for (int i = 0; i < 8; i++)
        ofs6 << hash_seed_1[i] + minminPosition_ << std::endl;
    ofs6.close();

    std::ofstream ofs7("c:\\temp\\LRc_quad2TMP.txt");
    for (int i = 0; i < 8; i++)
        ofs7 << hash_seed_2[i] + minminPosition_ << std::endl;
    ofs7.close();

#endif

    return true;
}

// ============================================================================================
// generateRefinedModel: compute model using non-minimal set of samples
// default operation is to use a weight of 1 for every data point
// ============================================================================================
bool Rigid3dEstimator::generateRefinedModel(std::vector<unsigned int>& sample,
                                          unsigned int numPoints,
                                          bool weighted,
                                          double* weights)
{
    // copy all 3D point pairs into eign
    opengv::points_t pts1, pts2;
    for (int i = 0; i < numPoints; i++)
    {
        opengv::point_t p1, p2;
        getPointPair3D(sample[i], p1, p2);
        pts1.push_back(p1);
        pts2.push_back(p2);
    }

    // solving using arun's algorithm
    opengv::transformation_t solution;
    solution = opengv::math::arun_complete(pts1, pts2);
    rotation_ = solution.block<3, 3>(0, 0);
    translation_ = solution.col(3);

    memcpy((void*)models_denorm_[0], solution.data(), model_size_ * sizeof(double));

    return true;
}


// ============================================================================================
// validateSample: check if minimal sample is valid
// ============================================================================================
bool Rigid3dEstimator::validateSample()
{
    // check (approximate) triangle congruency
    for (int i = 0; i < 3; i++)
    {
        double* pt1 = input_points_denorm_ + 6 * min_sample_[i];
        double* pt2 = input_points_denorm_ + 6 * min_sample_[(i+1)%3];

        double dist1 = DistR3(pt1, pt2);
        double dist2 = DistR3(pt1+3, pt2+3);
        if ((dist1 / dist2) < congruency_threshold_)
            return false;
        if ((dist2 / dist1) < congruency_threshold_)
            return false;
    }

   return true;
}


// ============================================================================================
// validateModel: check if model computed from minimal sample is valid
// ============================================================================================
bool Rigid3dEstimator::validateModel(const unsigned int modelIndex)
{
    return true;
}

double Rigid3dEstimator::getResidualError(double* correspondence, double* R_prime, double* T)
{
    double err, R_p2[3], R_p2_T[3], res;
    double* pt1 = correspondence;
    double* pt2 = correspondence + 3;

    MathTools::leftvmul(R_p2, pt2, R_prime, 3);
    MathTools::addvecs(R_p2_T, R_p2, T, 3);

    return DistR3(R_p2_T, pt1);
}

inline bool Rigid3dEstimator::checkSampsonVsThreshold(double* model, unsigned int* samplesInds, double threshold)
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
bool Rigid3dEstimator::evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested)
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
void Rigid3dEstimator::testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel)
{
    *degenerateModel = false;
    *upgradeModel = false;
}

// ============================================================================================
// upgradeDegenerateModel: try to upgrade degenerate model to non-degenerate by sampling from
// the set of outliers to the degenerate model
// ============================================================================================
unsigned int Rigid3dEstimator::upgradeDegenerateModel()
{
    return 0;
}


// ============================================================================================
// findWeights: given model and points, compute weights to be used in local optimization
// ============================================================================================
void Rigid3dEstimator::findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers,
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
void Rigid3dEstimator::storeModel(const unsigned int modelIndex, unsigned int numInliers)
{
    // save the current model as the best solution so far
    for (unsigned int i = 0; i < model_size_; ++i)
    {
        final_model_params_[i] = *(models_denorm_[modelIndex]+i);
    }
}

#endif
