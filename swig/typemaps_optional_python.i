/******************************************************************************
 *
 * SWIG typemaps for std::optional<T> (Python target)
 *
 * Handles automatic conversion between Python objects (or None) and
 * std::optional<T> types in C++.
 *
 *****************************************************************************/

%{
#include <optional>
%}

/**
 * @brief Macro helper to factorize SWIG typemaps for simple numeric integral types.
 * Converts Py_None to std::nullopt and valid Python integers to std::optional<TYPE>.
 */
%define OPTIONAL_TYPEMAP_NUMERIC(TYPE, CHECK_FUNC, CONV_FUNC, TYPE_NAME)
%typemap(in) std::optional<TYPE>
{
    if ($input == Py_None)
    {
        $1 = std::nullopt;
    }
    else if (CHECK_FUNC($input))
    {
        $1 = static_cast<TYPE>(CONV_FUNC($input));
    }
    else
    {
        SWIG_exception_fail(
            SWIG_TypeError,
            "Expected " TYPE_NAME " or None for optional<" TYPE_NAME ">"
        );
    }
}
%enddef

OPTIONAL_TYPEMAP_NUMERIC(int, PyLong_Check, PyLong_AsLong, "int")
OPTIONAL_TYPEMAP_NUMERIC(long, PyLong_Check, PyLong_AsLong, "int")
OPTIONAL_TYPEMAP_NUMERIC(long long, PyLong_Check, PyLong_AsLongLong, "int")

/**
 * @brief Macro helper to factorize SWIG typemaps for floating point types.
 * Accepts both Python floats and Python integers, converting Py_None to std::nullopt.
 */
%define OPTIONAL_TYPEMAP_FLOAT(TYPE, CONV_FUNC, TYPE_NAME)
%typemap(in) std::optional<TYPE>
{
    if ($input == Py_None)
    {
        $1 = std::nullopt;
    }
    else if (PyFloat_Check($input) || PyLong_Check($input))
    {
        $1 = static_cast<TYPE>(CONV_FUNC($input));
    }
    else
    {
        SWIG_exception_fail(
            SWIG_TypeError,
            "Expected float or None for optional<" TYPE_NAME ">"
        );
    }
}
%enddef

OPTIONAL_TYPEMAP_FLOAT(double, PyFloat_AsDouble, "double")
OPTIONAL_TYPEMAP_FLOAT(float, PyFloat_AsDouble, "float")

/**
 * @brief Typemap for std::optional<bool>.
 * Converts Python None to std::nullopt and Python booleans to std::optional<bool>.
 */
%typemap(in) std::optional<bool>
{
    if ($input == Py_None)
    {
        $1 = std::nullopt;
    }
    else if (PyBool_Check($input))
    {
        $1 = (PyObject_IsTrue($input) != 0);
    }
    else
    {
        SWIG_exception_fail(
            SWIG_TypeError,
            "Expected bool or None for optional<bool>"
        );
    }
}

/**
 * @brief Typemap for std::optional<std::string>.
 * Converts Python None to std::nullopt and UTF-8 Python strings to std::optional<std::string>.
 */
%typemap(in) std::optional<std::string>
{
    if ($input == Py_None)
    {
        $1 = std::nullopt;
    }
    else if (PyUnicode_Check($input))
    {
        Py_ssize_t size = 0;
        const char* str = PyUnicode_AsUTF8AndSize($input, &size);

        if (str == nullptr)
        {
            SWIG_exception_fail(
                SWIG_RuntimeError,
                "Cannot convert Python string"
            );
        }

        $1 = std::string(str, static_cast<size_t>(size));
    }
    else
    {
        SWIG_exception_fail(
            SWIG_TypeError,
            "Expected string or None for optional<string>"
        );
    }
}
