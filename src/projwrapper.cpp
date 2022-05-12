/*
 *
 */

#include <../include/projwrapper.hpp>
#include <iostream>

using projplot::PJContainer;
using projplot::ProjWrapper;
using projplot::ProjError;
using projplot::xy_t;

/**************************
 *
 *    ProjError
 *
 **************************/

ProjError::ProjError(const char* msg) : msg(msg)
{}

const char* ProjError::what() const noexcept
{
	return msg;
}




/**************************
 *
 *    PJContainer.
 *
 **************************/


PJContainer::PJContainer(const char* proj_str)
   : context(nullptr), projection(nullptr), _a(0.0), _f(0.0)
{
	std::cout << "proj_str: '" << proj_str << "'\n" << std::flush;
	/* Create the multi-threading context: */
	context = proj_context_create();
	if (!context){
		throw ProjError("Could not create multithreading context.");
	}

	/* Create the projection: */
	projection = proj_create(context, proj_str);
	if (!projection){
		throw ProjError("Could not create PROJ string.");
	}

	/* Get the ellipsoid through the CRS: */
	PJ* crs = proj_get_source_crs(context, projection);
	if (!crs){
		/* Assume GRS80 like PROJ default. */
		_a = 6378137.0;
		_f = 1.0 / 298.257222101;
	} else {
		PJ* ellps = proj_get_ellipsoid(context, crs);
		if (!ellps){
			proj_destroy(crs);
			throw ProjError("Could not obtain ellipsoid.");
		}

		/* Get the ellipsoid parameters: */
		double rf;
		int rc = proj_ellipsoid_get_parameters(context, ellps, &_a, NULL, NULL,
			                                   &rf);

		/* Cleanup: */
		proj_destroy(ellps);
		proj_destroy(crs);

		/* Check success: */
		if (!rc){
			throw ProjError("Could not obtain ellipsoid parameters.");
		}

		_f = 1.0 / rf;
	}

}


PJContainer::~PJContainer() {
	if (projection){
		proj_destroy(projection);
	}
	if (context){
		proj_context_destroy(context);
	}
}


const PJ* PJContainer::get() const
{
	return projection;
}


double PJContainer::a() const
{
	return _a;
}


double PJContainer::f() const
{
	return _f;
}





ProjWrapper::ProjWrapper(const char* proj_str)
   : proj_source(std::make_shared<PJContainer>(proj_str)),
     context(nullptr), workhorse(nullptr)
{
	/* Create the multi-threading context: */
	context = proj_context_create();
	if (!context){
		throw std::runtime_error("Could not create multithreading context.");
	}

	/* Create the projection: */
	workhorse = proj_clone(context, proj_source->get());
	if (!workhorse){
		throw std::runtime_error("Could not create PROJ string.");
	}
}


ProjWrapper::ProjWrapper(const ProjWrapper& other)
   : proj_source(other.proj_source), context(nullptr), workhorse(nullptr)
{
	/* Early sanity: */
	if (!proj_source->get())
		throw std::runtime_error("Trying to clone a wrapper with invalid "
		                         "context or projection");

	/* Create the multi-threading context: */
	context = proj_context_create();
	if (!context){
		throw std::runtime_error("Could not create multithreading context.");
	}

	/* Clone the projection: */
	workhorse = proj_clone(context, proj_source->get());
	if (!workhorse){
		throw std::runtime_error("Could not create PROJ string.");
	}
}


ProjWrapper::~ProjWrapper() {
	if (workhorse){
		proj_destroy(workhorse);
	}
	if (context){
		proj_context_destroy(context);
	}
}


void ProjWrapper::project(double lon, double lat, double& x, double& y) const
{
	PJ_COORD from(proj_coord(lon, lat, 0.0, 0.0));
	PJ_COORD to(proj_trans(workhorse, PJ_FWD, from));
	x = to.xy.x;
	y = to.xy.y;
}

xy_t ProjWrapper::project(double lon, double lat) const
{
	PJ_COORD from(proj_coord(lon, lat, 0.0, 0.0));
	PJ_COORD to(proj_trans(workhorse, PJ_FWD, from));
	return {to.xy.x, to.xy.y};
	//return {lon, lat};
}

xy_t ProjWrapper::project(const geo_t& lola) const
{
	PJ_COORD from(proj_coord(lola.lambda, lola.phi, 0.0, 0.0));
	PJ_COORD to(proj_trans(workhorse, PJ_FWD, from));
	return {to.xy.x, to.xy.y};
}


double ProjWrapper::a() const
{
	return proj_source->a();
}


double ProjWrapper::f() const
{
	return proj_source->f();
}
