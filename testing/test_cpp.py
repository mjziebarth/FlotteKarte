# Test the C++ code.
#
# Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
#
# Copyright (C) 2024 Technische Universität München
#
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
# the European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Licence is distributed on an "AS IS" basis,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and
# limitations under the Licence.

import os
from pathlib import Path
from subprocess import run


def test_compile():
    # Source directory:
    src_dir = Path(__file__).parent / 'cpp'

    # Ensure that the subprojects folder is linked:
    if not (src_dir / "subprojects").is_symlink():
        os.symlink(
            src_dir.parent.parent / "flottekarte" / "extensions"
                / "subprojects",
            src_dir / "subprojects"
        )

    # Setup the build dir:
    os.chdir(src_dir)
    builddir = src_dir / "builddir"
    if builddir.is_dir():
        run(("meson","setup","--reconfigure","builddir"), check=True)
    else:
        run(("meson","setup","builddir"), check=True)

    # Change to the builddir:
    os.chdir(builddir)

    # Run the compilation:
    run(("meson","compile"), check=True)

    # Run the executable:
    run((builddir / "test_flottekarte",), check=True)
