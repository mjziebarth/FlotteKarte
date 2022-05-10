/*
 *
 */

#include <../include/projwrapper.hpp>
#include <iostream>

using projplot::PJContainer;
using projplot::ProjWrapper;
using projplot::xy_t;

/**************************
 *
 *    PJContainer.
 *
 **************************/


PJContainer::PJContainer(const char* proj_str)
   : context(nullptr), projection(nullptr)
{
	std::cout << "proj_str: '" << proj_str << "'\n" << std::flush;
	/* Create the multi-threading context: */
	context = proj_context_create();
	if (!context){
		throw std::runtime_error("Could not create multithreading context.");
	}

	/* Create the projection: */
	projection = proj_create(context, proj_str);
	if (!projection){
		throw std::runtime_error("Could not create PROJ string.");
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
