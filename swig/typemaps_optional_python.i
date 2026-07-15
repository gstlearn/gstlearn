/******************************************************************************
 *
 * SWIG typemaps for std::optional<T>
 *
 *****************************************************************************/

%{
#include <optional>
%}


/*
 * optional<double>
 */
%typemap(in) std::optional<double>
{
    if ($input == Py_None)
    {
        $1 = std::nullopt;
    }
    else if (PyFloat_Check($input) || PyLong_Check($input))
    {
        $1 = static_cast<double>(PyFloat_AsDouble($input));
    }
    else
    {
        SWIG_exception_fail(
            SWIG_TypeError,
            "Expected float or None for optional<double>"
        );
    }
}


/*
 * optional<float>
 */
%typemap(in) std::optional<float>
{
    if ($input == Py_None)
    {
        $1 = std::nullopt;
    }
    else if (PyFloat_Check($input) || PyLong_Check($input))
    {
        $1 = static_cast<float>(PyFloat_AsDouble($input));
    }
    else
    {
        SWIG_exception_fail(
            SWIG_TypeError,
            "Expected float or None for optional<float>"
        );
    }
}


/*
 * optional<int>
 */
%typemap(in) std::optional<int>
{
    if ($input == Py_None)
    {
        $1 = std::nullopt;
    }
    else if (PyLong_Check($input))
    {
        $1 = static_cast<int>(PyLong_AsLong($input));
    }
    else
    {
        SWIG_exception_fail(
            SWIG_TypeError,
            "Expected int or None for optional<int>"
        );
    }
}


/*
 * optional<long>
 */
%typemap(in) std::optional<long>
{
    if ($input == Py_None)
    {
        $1 = std::nullopt;
    }
    else if (PyLong_Check($input))
    {
        $1 = static_cast<long>(PyLong_AsLong($input));
    }
    else
    {
        SWIG_exception_fail(
            SWIG_TypeError,
            "Expected int or None for optional<long>"
        );
    }
}


/*
 * optional<bool>
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


/*
 * optional<std::string>
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
