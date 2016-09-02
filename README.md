This is minimal Python bindings to C++ version of nvector library.

See

https://en.wikipedia.org/wiki/N-vector
http://www.navlab.net/nvector/
https://pypi.python.org/pypi/nvector

The latter was too bloated and dependency-ridden for my taste, so I did this.


Depends on: n_vector_library (C++ one at the nvector page) and eigen3.

n_vector_library is included here (it's 2-clause BSD I think).

hg clone eigen3 lib into eigen3 directory. It's header only, so no need to do anything else.

My code is 2-clause BSD or something like that.

