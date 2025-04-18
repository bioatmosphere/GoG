project('test','fortran',
    version:'0.1')

# Import the Python module
py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

#add_project_arguments('-arch arm64', language: 'fortran')

fortran_sources = files('Constants.f90')


# Run f2py to generate the wrapper
wrapper = custom_target('generate_wrapper',
    input: fortran_sources,
    output: 'constants.pyf',
    command: [
        'f2py',
        '--overwrite-signature',
        '-m', 'constants',
        '-h', '@OUTPUT@',
        '@INPUT@']
)

# Build the Python extension module using the generated wrapper
py_extension = py.extension_module(
    'constants',
    sources: [wrapper, fortran_sources],
    dependencies : py_dep,
    install: false,
    install_dir: py.get_install_dir()
)

# Print a success message
message('Fortran wrapper successfully built using f2py and Meson!')