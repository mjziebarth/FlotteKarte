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

#include <../include/types.hpp>
#include <../include/griddedinverter.hpp>
#include <memory>

#ifndef FLOTTEKARTE_AUGMENTEDPROJ_HPP
#define FLOTTEKARTE_AUGMENTEDPROJ_HPP

namespace flottekarte {

class AugmentedProj {
public:
	AugmentedProj(const char* proj_str);

	AugmentedProj(const AugmentedProj& other) = default;
	AugmentedProj(AugmentedProj&& other) = default;

	void project(double lam, double phi, double& x, double& y) const;
	xy_t project(double lam, double phi) const;
	xy_t project(const geo_t& lola) const;
	xy_t project(const geo_degrees_t& lola) const;

	geo_t inverse(const xy_t& xy_t) const;
	void inverse(double x, double y, double& lam, double& phi) const;

	double a() const;
	double f() const;

	const ProjWrapper& wrapper() const;

private:
	const ProjWrapper proj;
	const bool has_inverse;
	mutable std::shared_ptr<GriddedInverter> inverter;

	geo_t inverse_gd(const xy_t& xy) const;

};

}

#endif