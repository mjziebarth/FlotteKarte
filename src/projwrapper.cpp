/*
 *
 */

#include <../include/projwrapper.hpp>
#include <iostream>
#include <string>
#include <cstring>

using flottekarte::PJContainer;
using flottekarte::ProjWrapper;
using flottekarte::ProjError;
using flottekarte::xy_t;
using flottekarte::geo_t;

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

/* A very simple parser to obtain information about a and f from
 * a proj string. */
static void parse_proj_str_for_a_f(const char* proj_str, double& a, double& f)
{
	if (!proj_str)
		return;

	const char* c = proj_str;
	// 10000 characters should be fairly much enough for
	// a proj string.
	const size_t maxlen = strnlen(c, 10000);

	if (maxlen < 2){
		/* No possibility to place a keyword here. */
		return;
	}

	auto readkeyword = [&](char key, double& to) {
		size_t i0=0;
		for (size_t i=0; i<maxlen-2; ++i){
			if (c[i] == '\0')
				break;
			if ((i == 0) && (c[i] == key) && (c[i+1] == '=')){
				i0 = 2;
				break;
			} else if ((c[i] == '+' || c[i] == ' ') && c[i+1] == key
			           && c[i+2] == '=')
			{
				i0 = i+3;
				break;
			}
		}
		size_t i1 = maxlen; // maxlen >= 2 as per check above
		for (size_t i=i0; i<maxlen; ++i){
			if (c[i] == '\0' || c[i] == ' '){
				i1 = i;
				break;
			}
		}
		if ((i1 - i0) > 0) {
			try {
				std::string s(c+i0, i1-i0);
				double d = std::stod(s);
				if (d > 0)
					to = d;
			} catch (...) {
				/* Do nothing. */
			}
		}
	};

	readkeyword('a',a);
	readkeyword('f',f);
}



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

		/* Finally, try reading the ellipsoid parameters directly.
		 * If they are given as a=.. and f=.. in a PROJ string,
		 * no source CRS will be created for the PJ instance.
		 * So, `projection` will know the ellipsoid parameters
		 * but there is (likely?) no standard API way to obtain
		 * them from PJ.
		 * Instead, parse the PROJ string directly:
		 */
		parse_proj_str_for_a_f(proj_str, _a, _f);

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

xy_t ProjWrapper::project(const geo_degrees_t& lola) const
{
	return project(lola.to_radians());
}

geo_t ProjWrapper::inverse(const xy_t& xy) const
{
	if (proj_errno(workhorse))
		throw ProjError("Error state set before inverse.");
	PJ_COORD from(proj_coord(xy.x, xy.y, 0.0, 0.0));
	PJ_COORD to(proj_trans(workhorse, PJ_INV, from));
	int err = proj_errno(workhorse);
	if (err){
		/* Throw error with indicating the PROJ error.
		 * The thread-safe proj_context_errno_string routine was added
		 * only in PROJ version 8.0.0. If the local version is older,
		 * use the unsafe variant.
		 */
		#if PROJ_VERSION_MAJOR >= 8
		throw ProjError(proj_context_errno_string(context, err));
		#else
		throw ProjError(proj_errno_string(err));
		#endif
	}
	return {to.lp.lam, to.lp.phi};
}


double ProjWrapper::a() const
{
	return proj_source->a();
}


double ProjWrapper::f() const
{
	return proj_source->f();
}
