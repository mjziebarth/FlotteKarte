/*
 * Augmented version of the ProjWrapCpp ProjWrapper class.
 *
 * Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
 *
 * Copyright (C) 2022 Malte J. Ziebarth
 *
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
 * the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 */

#include <../include/augmentedproj.hpp>
#include <../include/invert.hpp>

using flottekarte::xy_t;
using flottekarte::geo_t;
using flottekarte::ProjWrapper;
using flottekarte::geo_degrees_t;
using flottekarte::AugmentedProj;
using flottekarte::GriddedInverter;
using flottekarte::gradient_descent_inverse_project;

AugmentedProj::AugmentedProj(const char* proj_str)
   : ProjWrapper(proj_str), has_inverse(ProjWrapper::has_inverse())
{
	if (!has_inverse){
		inverter = std::make_shared<GriddedInverter>(*this, 100, 50);
		if (!inverter)
			throw std::runtime_error("Could not initialize gridded inverter.");
	}
}


geo_t AugmentedProj::inverse(const xy_t& xy) const
{
	if (!has_inverse)
		return inverse_gd(xy);
	try {
		return ProjWrapper::inverse(xy);
	} catch (const ProjError& err) {
		return inverse_gd(xy);
	}
}


void AugmentedProj::inverse(double x, double y, double& lam, double& phi) const
{
	if (!has_inverse){
		geo_t geo(inverse_gd(xy_t(x,y)));
		lam = geo.lambda;
		phi = geo.phi;
	} else {
		try {
			ProjWrapper::inverse(x,y,lam,phi);
		} catch (const ProjError& err) {
			geo_t geo(inverse_gd(xy_t(x,y)));
			lam = geo.lambda;
			phi = geo.phi;
		}
	}
}


geo_t AugmentedProj::inverse_gd(const xy_t& xy) const
{
	/* Make sure the inverter is initialized: */
	if (!inverter){
		inverter = std::make_shared<GriddedInverter>(*this, 100, 50);
		if (!inverter)
			throw std::runtime_error("Could not initialize gridded inverter.");
	}

	/* Compute starting value using the inverter: */
	geo_t geo((*inverter)(xy));

	/* Use gradient descent to finish the job: */
	return gradient_descent_inverse_project(*this, xy, geo.lambda, geo.phi);
}
