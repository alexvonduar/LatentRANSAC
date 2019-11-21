#include "PnpFunctions.h"

namespace PTools
{

    void p3p_kneip_RT(
            const opengv::bearingVectors_t & f,
            const opengv::points_t & p,
        opengv::rotations_t & rotations,
        opengv::translations_t & translations,
        opengv::transformations_t & solutions)
    {
        opengv::point_t P1 = p[0];
        opengv::point_t P2 = p[1];
        opengv::point_t P3 = p[2];

        Eigen::Vector3d temp1 = P2 - P1;
        Eigen::Vector3d temp2 = P3 - P1;

        if (temp1.cross(temp2).norm() == 0)
            return;

        opengv::bearingVector_t f1 = f[0];
        opengv::bearingVector_t f2 = f[1];
        opengv::bearingVector_t f3 = f[2];

        Eigen::Vector3d e1 = f1;
        Eigen::Vector3d e3 = f1.cross(f2);
        e3 = e3 / e3.norm();
        Eigen::Vector3d e2 = e3.cross(e1);

        opengv::rotation_t T;
        T.row(0) = e1.transpose();
        T.row(1) = e2.transpose();
        T.row(2) = e3.transpose();

        f3 = T*f3;

        if (f3(2, 0) > 0)
        {
            f1 = f[1];
            f2 = f[0];
            f3 = f[2];

            e1 = f1;
            e3 = f1.cross(f2);
            e3 = e3 / e3.norm();
            e2 = e3.cross(e1);

            T.row(0) = e1.transpose();
            T.row(1) = e2.transpose();
            T.row(2) = e3.transpose();

            f3 = T*f3;

            P1 = p[1];
            P2 = p[0];
            P3 = p[2];
        }

        Eigen::Vector3d n1 = P2 - P1;
        n1 = n1 / n1.norm();
        Eigen::Vector3d n3 = n1.cross(P3 - P1);
        n3 = n3 / n3.norm();
        Eigen::Vector3d n2 = n3.cross(n1);

        opengv::rotation_t N;
        N.row(0) = n1.transpose();
        N.row(1) = n2.transpose();
        N.row(2) = n3.transpose();

        P3 = N*(P3 - P1);

        double d_12 = temp1.norm();
        double f_1 = f3(0, 0) / f3(2, 0);
        double f_2 = f3(1, 0) / f3(2, 0);
        double p_1 = P3(0, 0);
        double p_2 = P3(1, 0);

        double cos_beta = f1.dot(f2);
        double b = 1 / (1 - pow(cos_beta, 2)) - 1;

        if (cos_beta < 0)
            b = -sqrt(b);
        else
            b = sqrt(b);

        double f_1_pw2 = pow(f_1, 2);
        double f_2_pw2 = pow(f_2, 2);
        double p_1_pw2 = pow(p_1, 2);
        double p_1_pw3 = p_1_pw2 * p_1;
        double p_1_pw4 = p_1_pw3 * p_1;
        double p_2_pw2 = pow(p_2, 2);
        double p_2_pw3 = p_2_pw2 * p_2;
        double p_2_pw4 = p_2_pw3 * p_2;
        double d_12_pw2 = pow(d_12, 2);
        double b_pw2 = pow(b, 2);

        Eigen::Matrix<double, 5, 1> factors;

        factors(0, 0) = -f_2_pw2*p_2_pw4
            - p_2_pw4*f_1_pw2
            - p_2_pw4;

        factors(1, 0) = 2 * p_2_pw3*d_12*b
            + 2 * f_2_pw2*p_2_pw3*d_12*b
            - 2 * f_2*p_2_pw3*f_1*d_12;

        factors(2, 0) = -f_2_pw2*p_2_pw2*p_1_pw2
            - f_2_pw2*p_2_pw2*d_12_pw2*b_pw2
            - f_2_pw2*p_2_pw2*d_12_pw2
            + f_2_pw2*p_2_pw4
            + p_2_pw4*f_1_pw2
            + 2 * p_1*p_2_pw2*d_12
            + 2 * f_1*f_2*p_1*p_2_pw2*d_12*b
            - p_2_pw2*p_1_pw2*f_1_pw2
            + 2 * p_1*p_2_pw2*f_2_pw2*d_12
            - p_2_pw2*d_12_pw2*b_pw2
            - 2 * p_1_pw2*p_2_pw2;

        factors(3, 0) = 2 * p_1_pw2*p_2*d_12*b
            + 2 * f_2*p_2_pw3*f_1*d_12
            - 2 * f_2_pw2*p_2_pw3*d_12*b
            - 2 * p_1*p_2*d_12_pw2*b;

        factors(4, 0) = -2 * f_2*p_2_pw2*f_1*p_1*d_12*b
            + f_2_pw2*p_2_pw2*d_12_pw2
            + 2 * p_1_pw3*d_12
            - p_1_pw2*d_12_pw2
            + f_2_pw2*p_2_pw2*p_1_pw2
            - p_1_pw4
            - 2 * f_2_pw2*p_2_pw2*p_1*d_12
            + p_2_pw2*f_1_pw2*p_1_pw2
            + f_2_pw2*p_2_pw2*d_12_pw2*b_pw2;

        std::vector<double> realRoots = opengv::math::o4_roots(factors);

        for (int i = 0; i < 4; i++)
        {
            double cot_alpha =
                (-f_1*p_1 / f_2 - realRoots[i] * p_2 + d_12*b) /
                (-f_1*realRoots[i] * p_2 / f_2 + p_1 - d_12);

            double cos_theta = realRoots[i];
            double sin_theta = sqrt(1 - pow(realRoots[i], 2));
            double sin_alpha = sqrt(1 / (pow(cot_alpha, 2) + 1));
            double cos_alpha = sqrt(1 - pow(sin_alpha, 2));

            if (cot_alpha < 0)
                cos_alpha = -cos_alpha;

            opengv::translation_t C;
            C(0, 0) = d_12*cos_alpha*(sin_alpha*b + cos_alpha);
            C(1, 0) = cos_theta*d_12*sin_alpha*(sin_alpha*b + cos_alpha);
            C(2, 0) = sin_theta*d_12*sin_alpha*(sin_alpha*b + cos_alpha);

            C = P1 + N.transpose()*C;

            opengv::rotation_t R;
            R(0, 0) = -cos_alpha;
            R(0, 1) = -sin_alpha*cos_theta;
            R(0, 2) = -sin_alpha*sin_theta;
            R(1, 0) = sin_alpha;
            R(1, 1) = -cos_alpha*cos_theta;
            R(1, 2) = -cos_alpha*sin_theta;
            R(2, 0) = 0.0;
            R(2, 1) = -sin_theta;
            R(2, 2) = cos_theta;

            R = N.transpose()*R.transpose()*T;

            opengv::transformation_t solution;
            solution.col(3) = C;
            solution.block<3, 3>(0, 0) = R;

            solutions.push_back(solution);
            rotations.push_back(R);
            translations.push_back(C);
        }
    }

    void epnp_RT(
        const opengv::bearingVectors_t & f,
        const opengv::points_t & p,
        opengv::rotation_t & rotation,
        opengv::translation_t & translation,
        opengv::transformation_t & solution)
    {

        opengv::absolute_pose::modules::Epnp PnP;
        PnP.set_maximum_number_of_correspondences(p.size());
        PnP.reset_correspondences();

        for (unsigned int i = 0; i < p.size(); i++)
        {
            opengv::point_t p_i = p[i]; // adapter.getPoint(indices[i]);
            opengv::bearingVector_t f_i = f[i];

            PnP.add_correspondence(p_i[0], p_i[1], p_i[2], f_i[0], f_i[1], f_i[2]);
        }

        double R_epnp[3][3], t_epnp[3];
        PnP.compute_pose(R_epnp, t_epnp);

        for (int r = 0; r < 3; r++)
        {
            for (int c = 0; c < 3; c++)
                rotation(r, c) = R_epnp[r][c];
        }

        translation[0] = t_epnp[0];
        translation[1] = t_epnp[1];
        translation[2] = t_epnp[2];

        //take inverse transformation
        rotation.transposeInPlace();
        translation = -rotation * translation;

        solution.col(3) = translation;
        solution.block<3, 3>(0, 0) = rotation;

    }
}