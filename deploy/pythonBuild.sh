#!/bin/bash

if [ $2 = "py35" ]; then
	export PATH=/opt/python/p35/bin:$PATH
fi

if [ $2 = "py36" ]; then
	export PATH=/opt/python/p36/bin:$PATH
fi

if [ $2 = "py37" ]; then
	export PATH=/opt/python/p37/bin:$PATH
fi

if [ $2 = "py38" ]; then
	export PATH=/opt/python/p38/bin:$PATH
fi

if [ $2 = "py39" ]; then
	export PATH=/opt/python/p39/bin:$PATH
fi

if [ $2 = "py310" ]; then
	export PATH=/opt/python/p310/bin:$PATH
fi

cmake -Bbuild  -DHDF5_ROOT:PATH=/opt/hdf5-1.12.1 -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_C_FLAGS="-DH5_USE_110_API"
cmake --build build --target python_build -- -j8
cd $1
python3 setup.py bdist_wheel --python-tag=$2 --plat-name=manylinux_2_17_x86_64



