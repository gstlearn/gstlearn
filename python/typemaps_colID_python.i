/***************************************************************************/
/*                                                                         */
/*  Python typemaps for gstlrn::ColID                                      */
/*                                                                         */
/*  Dual-binding support for both:                                         */
/*    1. const gstlrn::ColID&  (Preferred API facade in Db.hpp)            */
/*    2. const gstlrn::ColID&& (Legacy / Low-level engine in DbData.hpp)   */
/*                                                                         */
/*  This typemap transparently converts high-level Python arguments into   */
/*  a C++ ColID temporary instance.                                        */
/*                                                                         */
/*  Accepted Python inputs:                                                */
/*                                                                         */
/*    - "name"                  -> ColID("name")                           */
/*    - "name", version         -> ColID("name", version)                  */
/*    - icol                    -> ColID(icol)                             */
/*    - icol, version           -> ColID(icol, version)                    */
/*    - RoleID                  -> ColID(RoleID)                           */
/*    - RoleID, version         -> ColID(RoleID, version)                  */
/*    - ERole                   -> ColID(ERole, 0)                         */
/*    - ERole, version          -> ColID(ERole, 0), version)               */
/*                                                                         */
/***************************************************************************/

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
    const gstlrn::ColID&,
    gstlrn::ColID&&
{
    PyObject *obj = $input;

    if (PyTuple_Check(obj))
    {
        if (PyTuple_GET_SIZE(obj) == 2)
        {
            obj = PyTuple_GET_ITEM(obj, 0);
        }
        else
        {
            $1 = 0;
            obj = nullptr;
        }
    }

    if (obj != nullptr)
    {
        if (PyUnicode_Check(obj) || PyLong_Check(obj))
        {
            $1 = 1;
        }
        else
        {
            static swig_type_info *type_ColID  = SWIG_TypeQuery("gstlrn::ColID *");
            static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
            static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

            void *ptr = nullptr;
            if (type_ERole && SWIG_IsOK(SWIG_ConvertPtr(obj, &ptr, type_ERole, 0)) && ptr)
            {
                $1 = 1;
            }
            else if (type_RoleID && SWIG_IsOK(SWIG_ConvertPtr(obj, &ptr, type_RoleID, 0)) && ptr)
            {
                $1 = 1;
            }
            else if (type_ColID && SWIG_IsOK(SWIG_ConvertPtr(obj, &ptr, type_ColID, 0)) && ptr)
            {
                $1 = 1;
            }
            else
            {
                $1 = 0;
            }
        }
    }
}


%typemap(in)
    const gstlrn::ColID&,
    gstlrn::ColID&&
{
    PyObject *obj = $input;
    gstlrn::Id version = 0;

    if (PyTuple_Check(obj))
    {
        if (PyTuple_GET_SIZE(obj) != 2)
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "ColID tuple must contain exactly two elements");
        }

        PyObject *first  = PyTuple_GET_ITEM(obj, 0);
        PyObject *second = PyTuple_GET_ITEM(obj, 1);

        if (!PyLong_Check(second))
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "Second tuple element must be an integer version");
        }

        version = static_cast<gstlrn::Id>(PyLong_AsLong(second));
        obj = first;
    }

    if (PyUnicode_Check(obj))
    {
        Py_ssize_t size = 0;
        const char *str = PyUnicode_AsUTF8AndSize(obj, &size);

        if (str == nullptr)
        {
            SWIG_exception_fail(
                SWIG_RuntimeError,
                "Cannot convert Python string");
        }

        $1 = new gstlrn::ColID(std::string_view(str, size), version);
    }
    else if (PyLong_Check(obj))
    {
        gstlrn::Id icol = static_cast<gstlrn::Id>(PyLong_AsLong(obj));
        $1 = new gstlrn::ColID(icol, version);
    }
    else
    {
        static swig_type_info *type_ColID  = SWIG_TypeQuery("gstlrn::ColID *");
        static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
        static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

        gstlrn::ERole *role = nullptr;
        int res_erole = SWIG_ConvertPtr(
            obj,
            (void**)&role,
            type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole,
            0);

        if (SWIG_IsOK(res_erole) && role != nullptr)
        {
            $1 = new gstlrn::ColID(*role, version);
        }
        else
        {
            gstlrn::RoleID *roleID = nullptr;
            int res_roleid = SWIG_ConvertPtr(
                obj,
                (void**)&roleID,
                type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID,
                0);

            if (SWIG_IsOK(res_roleid) && roleID != nullptr)
            {
                $1 = new gstlrn::ColID(*roleID, version);
            }
            else
            {
                gstlrn::ColID *colid = nullptr;
                int res_colid = SWIG_ConvertPtr(
                    obj,
                    (void**)&colid,
                    type_ColID ? type_ColID : SWIGTYPE_p_gstlrn__ColID,
                    0);

                if (SWIG_IsOK(res_colid) && colid != nullptr)
                {
                    if (version != 0)
                    {
                        $1 = new gstlrn::ColID(gstlrn::ColID::create(*colid, version));
                    }
                    else
                    {
                        $1 = new gstlrn::ColID(*colid);
                    }
                }
                else
                {
                    SWIG_exception_fail(
                        SWIG_TypeError,
                        "Expected string, integer, RoleID, ERole, "
                        "(string,version), (integer,version), "
                        "(RoleID,version) or (ERole,version)");
                }
            }
        }
    }
}


%typemap(freearg)
    const gstlrn::ColID&,
    gstlrn::ColID&&
{
    delete $1;
}
