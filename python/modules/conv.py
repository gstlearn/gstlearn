import gstlearn as gl
import numpy as np

try:
    import pandas as pd
    import scipy.sparse as sc
except ModuleNotFoundError as ex:
    msg = (
        "Python dependencies 'pandas' and 'scipy' not found.\n"
        "To install them alongside gstlearn, please run `pip install gstlearn[conv]'"
    )
    raise ModuleNotFoundError(msg) from ex


def Db_toTL(self, flagLocate=False):
    df = pd.DataFrame(
        self.getAllColumns().reshape(-1, self.getNSample()).T,
        columns=self.getAllNames(),
    )

    if flagLocate:
        for j, i in enumerate(self.getAllNames()):
            df[i].locator = self.getLocators()[j]
    return df


# TODO : This below (and all other setattr for toTL) overrides
# DECLARE_TOTL usage (not needed in python ?)
setattr(gl.Db, "toTL", Db_toTL)


def Db_fromPandas(df):
    # Create an empty Db
    dat = gl.Db()
    # And import all columns in one a loop using [] operator
    for field in df.columns:
        mycol = df[field]
        if mycol.dtype == "float64" or mycol.dtype == "int64":
            dat[field] = mycol
    return dat


gl.Db.fromTL = staticmethod(Db_fromPandas)


def table_toTL(self):
    # As a Panda Data Frame
    colnames = self.getColumnNames()
    rownames = self.getRowNames()
    if len(colnames) == 0:
        colnames = None
    if len(rownames) == 0:
        rownames = None
    Anp = pd.DataFrame(
        self.getValues(False).reshape(self.getNRows(), self.getNCols()),
        columns=colnames,
        index=rownames,
    )
    return Anp


setattr(gl.Table, "toTL", table_toTL)


def vario_toTL(self, idir, ivar, jvar):
    sw = self.getSwVec(idir, ivar, jvar, False)
    hh = self.getHhVec(idir, ivar, jvar, False)
    gg = self.getGgVec(idir, ivar, jvar, False, False, False)
    array = np.vstack((sw, hh, gg)).T
    colnames = np.array(["sw", "hh", "gg"])
    return pd.DataFrame(array, columns=colnames)


setattr(gl.Vario, "toTL", vario_toTL)


def vario_updateFromPanda(self, pf, idir, ivar, jvar):
    vario = self
    ndir = vario.getNDir()
    nvar = vario.getNVar()
    if idir < 0 or idir >= ndir:
        return vario
    if ivar < 0 or ivar >= nvar:
        return vario
    if jvar < 0 or jvar >= nvar:
        return vario
    nlag = vario.getNLagTotal(idir)
    if len(pf.index) != nlag:
        return vario

    vario.setSwVec(idir, ivar, jvar, pf["sw"])
    vario.setHhVec(idir, ivar, jvar, pf["hh"])
    vario.setGgVec(idir, ivar, jvar, pf["gg"])
    return vario


setattr(gl.Vario, "updateFromPanda", vario_updateFromPanda)


def matrix_toTL(self):
    if self.isSparse():
        NF_T = self.getMatrixToTriplet()
        return Triplet_toTL(NF_T)
    else:
        return np.array(self.getValues(False)).reshape(self.getNRows(), self.getNCols())


setattr(gl.MatrixDense, "toTL", matrix_toTL)
setattr(gl.MatrixSquare, "toTL", matrix_toTL)
setattr(gl.MatrixSymmetric, "toTL", matrix_toTL)
setattr(gl.MatrixSparse, "toTL", matrix_toTL)
setattr(gl.ProjMatrix, "toTL", matrix_toTL)
setattr(gl.PrecisionOpMultiMatrix, "toTL", matrix_toTL)
setattr(gl.ProjMultiMatrix, "toTL", matrix_toTL)


def Triplet_toTL(self):
    return sc.csc_matrix(
        (
            np.array(self.getValues()),
            (np.array(self.getRows()), np.array(self.getCols())),
        ),
        shape=(self.getNRows() + 1, self.getNCols() + 1),
    )


setattr(gl.NF_Triplet, "toTL", Triplet_toTL)
