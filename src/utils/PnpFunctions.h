#ifndef PTOOLS_H
#define PTOOLS_H

#include "MathFunctions.h"
#include <opengv/types.hpp>
#include <opengv/absolute_pose/modules/main.hpp>
#include <opengv/absolute_pose/modules/Epnp.hpp>
#include <opengv/math/roots.hpp>

namespace PTools
{
	void p3p_kneip_RT(
		const opengv::bearingVectors_t & f,
		const opengv::points_t & p,
		opengv::rotations_t & rotations,
		opengv::translations_t & translations,
		opengv::transformations_t & solutions);

	void epnp_RT(
		const opengv::bearingVectors_t & f,
		const opengv::points_t & p,
		opengv::rotation_t & rotation,
		opengv::translation_t & translation,
		opengv::transformation_t & solution);
}
#endif