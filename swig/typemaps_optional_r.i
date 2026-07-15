/******************************************************************************
 *
 * SWIG typemaps for std::optional<T>
 *
 *****************************************************************************/
/******************************************************************************
 *
 * SWIG R typemaps for std::optional<T>
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
    if (Rf_isNull($input))
    {
        $1 = std::nullopt;
    }
    else if (Rf_isReal($input) || Rf_isInteger($input))
    {
        $1 = static_cast<double>(Rf_asReal($input));
    }
    else
    {
        Rf_error(
            "Expected numeric or NULL for optional<double>");
    }
}


/*
 * optional<float>
 */
%typemap(in) std::optional<float>
{
    if (Rf_isNull($input))
    {
        $1 = std::nullopt;
    }
    else if (Rf_isReal($input) || Rf_isInteger($input))
    {
        $1 = static_cast<float>(Rf_asReal($input));
    }
    else
    {
        Rf_error(
            "Expected numeric or NULL for optional<float>");
    }
}


/*
 * optional<int>
 */
%typemap(in) std::optional<int>
{
    if (Rf_isNull($input))
    {
        $1 = std::nullopt;
    }
    else if (Rf_isInteger($input))
    {
        $1 = INTEGER($input)[0];
    }
    else
    {
        Rf_error(
            "Expected integer or NULL for optional<int>");
    }
}


/*
 * optional<long>
 */
%typemap(in) std::optional<long>
{
    if (Rf_isNull($input))
    {
        $1 = std::nullopt;
    }
    else if (Rf_isInteger($input))
    {
        $1 = static_cast<long>(INTEGER($input)[0]);
    }
    else
    {
        Rf_error(
            "Expected integer or NULL for optional<long>");
    }
}


/*
 * optional<bool>
 */
%typemap(in) std::optional<bool>
{
    if (Rf_isNull($input))
    {
        $1 = std::nullopt;
    }
    else if (Rf_isLogical($input))
    {
        $1 = (LOGICAL($input)[0] != 0);
    }
    else
    {
        Rf_error(
            "Expected logical or NULL for optional<bool>");
    }
}


/*
 * optional<std::string>
 */
%typemap(in) std::optional<std::string>
{
    if (Rf_isNull($input))
    {
        $1 = std::nullopt;
    }
    else if (Rf_isString($input))
    {
        $1 = std::string(
            CHAR(STRING_ELT($input,0)));
    }
    else
    {
        Rf_error(
            "Expected character or NULL for optional<string>");
    }
}
