
#ifndef GSTLEARN_EXPORT_H
#define GSTLEARN_EXPORT_H

#ifdef GSTLEARN_STATIC_DEFINE
#  define GSTLEARN_EXPORT
#  define GSTLEARN_NO_EXPORT
#else
#  ifndef GSTLEARN_EXPORT
#    ifdef shared_EXPORTS
        /* We are building this library */
#      define GSTLEARN_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define GSTLEARN_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef GSTLEARN_NO_EXPORT
#    define GSTLEARN_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef GSTLEARN_DEPRECATED
#  define GSTLEARN_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef GSTLEARN_DEPRECATED_EXPORT
#  define GSTLEARN_DEPRECATED_EXPORT GSTLEARN_EXPORT GSTLEARN_DEPRECATED
#endif

#ifndef GSTLEARN_DEPRECATED_NO_EXPORT
#  define GSTLEARN_DEPRECATED_NO_EXPORT GSTLEARN_NO_EXPORT GSTLEARN_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef GSTLEARN_NO_DEPRECATED
#    define GSTLEARN_NO_DEPRECATED
#  endif
#endif

#ifdef SWIG
#  undef GSTLEARN_EXPORT
#  undef GSTLEARN_NO_EXPORT
#  define GSTLEARN_EXPORT
#  define GSTLEARN_NO_EXPORT
#endif

#endif /* GSTLEARN_EXPORT_H */
