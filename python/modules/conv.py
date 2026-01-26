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
        # for j, i in enumerate(self.getAllNames()):
        #    df[i].locator = self.getLocators()[j]
        # This no longer works with pandas >= 3.0
        # because df[i] is a copy and not a view
        print("Warning: 'flagLocate' argument no longer supported")
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


def matrix_general_toLatex(self, col_titles=None, row_titles=None, precision=3):
    """
    Retourne une chaîne LaTeX compatible Jupyter Notebook pour une matrice avec ou sans titres.

    self : the input matrix (N x P)
    col_titles : liste de titres de colonnes (longueur P) ou None
    row_titles : liste de titres de lignes (longueur N) ou None
    precision : nombre de chiffres après la virgule
    """
    values = self.getValues()
    N = self.getNRows()
    P = self.getNCols()
    matrix = [values[i::P] for i in range(P)]

    # Déterminer l'alignement des colonnes
    # Si on a row_titles, ajouter un 'c' pour la première colonne
    # Déterminer l'alignement des colonnes
    num_cols = P + (1 if row_titles is not None else 0)
    col_format = "c" * num_cols  # alignement centré

    lines = []

    # Ligne d'en-tête si col_titles
    if col_titles is not None:
        if row_titles is not None:
            header = tuple([""]) + col_titles
        else:
            header = col_titles
        lines.append(" & ".join(header))

    # Lignes du tableau
    for i, row in enumerate(matrix):
        formatted_row = [f"{x:.{precision}f}" for x in row]
        if row_titles is not None:
            line = [row_titles[i]] + formatted_row
        else:
            line = formatted_row
        lines.append(" & ".join(line))

    # Concatenation avec \\ entre les lignes
    latex_body = " \\\\\n".join(lines)

    # Code final avec retours à la ligne
    latex_code = (
        "$$\\left[\\begin{array}{"
        + col_format
        + "}\n"
        + latex_body
        + "\n\\end{array}\\right]$$"
    )

    return latex_code


def matrix_toLatex(self, precision=3):
    return matrix_general_toLatex(self, None, None)


def table_toLatex(self, precision=3):
    colnames = self.getColumnNames()
    rownames = self.getRowNames()
    return matrix_general_toLatex(self, colnames, rownames, precision)


setattr(gl.MatrixDense, "toLatex", matrix_toLatex)
setattr(gl.MatrixSquare, "toLatex", matrix_toLatex)
setattr(gl.MatrixSymmetric, "toLatex", matrix_toLatex)
setattr(gl.MatrixSparse, "toLatex", matrix_toLatex)
setattr(gl.Table, "toLatex", table_toLatex)

setattr(gl.ProjMatrix, "toLatex", matrix_toLatex)
setattr(gl.PrecisionOpMultiMatrix, "toLatex", matrix_toLatex)
setattr(gl.ProjMultiMatrix, "toLatex", matrix_toLatex)
