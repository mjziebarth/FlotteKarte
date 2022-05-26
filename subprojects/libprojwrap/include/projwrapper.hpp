/*
 * Inverse projection based on inverting PROJ forward projection.
 */

#include <proj.h>
#include <memory>
#include <stdexcept>
#include <projwrappertypes.hpp>

#ifndef FLOTTEKARTE_PROJWRAPPER_HPP
#define FLOTTEKARTE_PROJWRAPPER_HPP

namespace projwrapper {

/* An exception class for PROJ errors: */
class ProjError : public std::exception
{
public:
	ProjError(const char* msg);

	virtual const char* what() const noexcept;

private:
	const char* msg;
};



/* A container for an originally created PJ object: */
class PJContainer {
public:
    PJContainer(const char* proj_str);
    ~PJContainer();
    const PJ* get() const;

	double a() const;
	double f() const;

private:
    PJ_CONTEXT* context;
    PJ* projection;
    double _a;
    double _f;
};


/* The main wrapper class. */
class ProjWrapper {
public:
    ProjWrapper(const char* proj_str);
    ProjWrapper(const ProjWrapper& other);
    ~ProjWrapper();

    void project(double lon, double lat, double& x, double& y) const;
    xy_t project(double lon, double lat) const;
	xy_t project(const geo_t& lola) const;
	xy_t project(const geo_degrees_t& lola) const;

	geo_t inverse(const xy_t& xy_t) const;

	double a() const;
	double f() const;

private:
    std::shared_ptr<PJContainer> proj_source;

    PJ_CONTEXT* context;
    PJ* workhorse;
};


}

#endif
