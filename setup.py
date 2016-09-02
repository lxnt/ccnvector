from distutils.core import setup, Extension

ccnvector = Extension('ccnvector',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['./n_vector_library/include', '.'],
                    libraries = ['m', 'stdc++'],
                    library_dirs = [],
                    sources = ['ccnvector.cc', './n_vector_library/src/n_vector_library.cpp'])

setup (name = 'ccnvector',
       version = '1.0',
       description = 'leanest iface to nvector',
       author = 'Stoddard',
       author_email = 'stoddard@lxnt.info',
       url = 'https://docs.python.org/extending/building',
       long_description = '''
...
''',
       ext_modules = [ccnvector])
