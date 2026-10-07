/***************************************************************************/
/*                                                                         */
/*  Python typemap for const gstlrn::RoleID&                               */
/*                                                                         */
/*  Accepted Python objects:                                               */
/*                                                                         */
/*      RoleID(...)             -> RoleID instance                         */
/*      (RoleID(...),index)     -> RoleID(RoleID(...), index)              */
/*                                                                         */
/*      ERole.X                 -> RoleID(ERole::X, 0)                     */
/*      (ERole.X,index)         -> RoleID(ERole::X, index)                 */
/*                                                                         */
/***************************************************************************/

%typemap(in) const gstlrn::RoleID& (gstlrn::RoleID temp)
{
    PyObject *obj = $input;
    gstlrn::Id index = 0;
    bool has_custom_index = false;

    /**********************************************************************/
    /* 1. Optional tuple: (object, index)                                 */
    /**********************************************************************/

    if (PyTuple_Check(obj))
    {
        if (PyTuple_GET_SIZE(obj) != 2)
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "RoleID tuple must contain exactly two elements: (RoleID|ERole, index)");
        }

        PyObject *first  = PyTuple_GET_ITEM(obj, 0);
        PyObject *second = PyTuple_GET_ITEM(obj, 1);

        if (!PyLong_Check(second))
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "Second tuple element must be an integer index");
        }

        index = static_cast<gstlrn::Id>(PyLong_AsLong(second));
        has_custom_index = true;
        obj = first;
    }

    /**********************************************************************/
    /* 2. C++ Wrapped Objects (RoleID or ERole)                           */
    /**********************************************************************/

    static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
    static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

    gstlrn::RoleID *roleID = nullptr;
    gstlrn::ERole  *role   = nullptr;

    // Check RoleID first
    int res1 = SWIG_ConvertPtr(
        obj,
        (void**)&roleID,
        type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID,
        0);

    if (SWIG_IsOK(res1) && roleID != nullptr)
    {
        temp = *roleID;
        if (has_custom_index)
        {
            temp.setIndex(index);
        }
        $1 = &temp;
    }
    else
    {
        // Fallback to ERole
        int res2 = SWIG_ConvertPtr(
            obj,
            (void**)&role,
            type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole,
            0);

        if (SWIG_IsOK(res2) && role != nullptr)
        {
            temp = gstlrn::RoleID(*role, index);
            $1 = &temp;
        }
        else
        {
            SWIG_exception_fail(
                SWIG_TypeError,
                "Expected RoleID, ERole, (RoleID, index) or (ERole, index)");
        }
    }
}

// Crucial step: Allow SWIG type-checking to recognize ERole and RoleID as valid inputs
%typecheck(SWIG_TYPECHECK_POINTER) const gstlrn::RoleID&
{
    static swig_type_info *type_RoleID = SWIG_TypeQuery("gstlrn::RoleID *");
    static swig_type_info *type_ERole  = SWIG_TypeQuery("gstlrn::ERole *");

    PyObject *obj = $input;
    if (PyTuple_Check(obj) && PyTuple_GET_SIZE(obj) == 2)
    {
        obj = PyTuple_GET_ITEM(obj, 0);
    }

     void *vptr = nullptr;
     int res1 = SWIG_ConvertPtr(obj, &vptr, type_RoleID ? type_RoleID : SWIGTYPE_p_gstlrn__RoleID, 0);
     int res2 = SWIG_ConvertPtr(obj, &vptr, type_ERole ? type_ERole : SWIGTYPE_p_gstlrn__ERole, 0);

     $1 = (SWIG_IsOK(res1) || SWIG_IsOK(res2)) ? 1 : 0;
}
