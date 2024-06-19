################################################################################
#                                                                              #
#                         gstlearn Python package                              #
#                                                                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://gstlearn.org                                                #
# License: BSD 3-clause                                                        #
#                                                                              #
################################################################################

import setuptools
from setuptools.command.build_ext import build_ext


class DummyExtensionBuild(build_ext):
    def run(self) -> None:
        return


setuptools.setup(
    # we add an empty dummy extension to get a platform wheel...
    ext_modules=[
        setuptools.Extension(
            name="dummy_extension_for_platform_wheel",
            sources=[],
        )
    ],
    # ...and we make sure to skip its compilation
    cmdclass={"build_ext": DummyExtensionBuild},
)
