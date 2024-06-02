/*
 * Test linear interpolator.
 *
 * Authors: Malte J. Ziebarth (malte.ziebarth@tum.de)
 *
 * Copyright (C) 2024 Technische Universität München
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

#include <interpolate.hpp>
#include <iostream>
#include <random>


static double fun(double x, double y)
{
    return 2.4 * x - 0.39 * y;
}

namespace flottekarte {


void test_interpolate()
{
    std::cout << "test_interpolate()\n" << std::flush;
    /* Generate the data: */
    size_t nx = 10;
    size_t ny = 8;
    std::vector<double> z(nx*ny, 0.0);

    double xmin = 0.3;
    double xmax = 2.8;
    double ymin = 2.1;
    double ymax = 3.4;

    std::cout << "setup array...\n" << std::flush;
    for (size_t i=0; i<nx; ++i){
        double x = xmin + ((xmax - xmin)*i) / (nx-1);
        for (size_t j=0; j<ny; ++j){
            double y = ymin + ((ymax - ymin)*j) / (ny-1);
            z[ny*i + j] = fun(x,y);
        }
    }

    std::cout << "setup interpolator...\n" << std::flush;

    /* Generate the interpolator: */
    LinearGrid2DInterpolator<false> interp(
        xmin, xmax, nx, ymin, ymax, ny, z.data()
    );

    std::cout << "test...\n" << std::flush;


    std::default_random_engine rng(8299832);
    std::uniform_real_distribution dist(0.0, 1.0);

    for (size_t i=0; i<30; ++i)
    {
        xy_t xy_i(
            xmin + (xmax - xmin) * dist(rng),
            ymin + (ymax - ymin) * dist(rng)
        );

        double z0 = fun(xy_i.x, xy_i.y);
        double z1 = interp(xy_i)[0];
        if (std::abs(z0-z1) > 1e-12){
            std::cerr << "[" << i << "]\n";
            std::cerr << "xy_i = (" << xy_i.x << ", " << xy_i.y << ")\n";
            std::cerr << "z0 = " << z0 << "\n"
                << "z1 = " << z1 << "\n";
            throw std::runtime_error("Results wrong.");
        }
    }

    std::cout << "done!\n" << std::flush;

}

} // end namespace