/*
 * Inverse projection based on inverting PROJ forward projection.
 */

#include <proj.h>
#include <memory>
#include <stdexcept>
#include <types.hpp>

#ifndef PROJPLOT_PROJWRAPPER_HPP
#define PROJPLOT_PROJWRAPPER_HPP

namespace projplot {

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

private:
    PJ_CONTEXT* context;
    PJ* projection;
};


/* The main wrapper class. */
class ProjWrapper {
public:
    ProjWrapper(const char* proj_str);
    ProjWrapper(const ProjWrapper& other);
    ~ProjWrapper();

    void project(double lon, double lat, double& x, double& y) const;
    xy_t project(double lon, double lat) const;

private:
    std::shared_ptr<PJContainer> proj_source;

    PJ_CONTEXT* context;
    PJ* workhorse;
};


}

#endif
