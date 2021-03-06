/* ----------------------------------------------------------------------------
    This header file includes class Param declaration.
    It reads and filters program arguments.
---------------------------------------------------------------------------- */

#ifndef PARAM_DEFINED
#define PARAM_DEFINED

class Param {
public:
    static const char *read(
            const int a_argc, const char **a_argv,  // arguments
            const char *a_param,    // param flag e.g. "-i"
            const char *a_def);     // default value if the flag is not found
};

#endif // PARAM_DEFINED

