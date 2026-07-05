import gstlearn as gl
import numpy as np
import pandas as pd

#
# This test is meant to check the manipulation of the new DataBase under Python
#

print("Instanciation de DbData")
data = gl.DbData()

print("Activation du test")
data.mytest()

print("Ajout d'une colonne de type Double")
data.addArrayDbl([1.0, 2.0, 3.0], "hello")
print("Ajout d'une colonne de type Int")
data.addArrayInt([5, 6, 7], "world")
print("Ajout d'une colonne de type String")
data.addArrayStr(["foo", "bar", "baz"], "foobar")
print("Presence de la colonne 'hello': ", data.hasArray("hello"))
print("Presence de la colonne 'world': ", data.hasArray("world"))
print("Presence de la colonne 'foobar': ", data.hasArray("foobar"))

print("Contenu de l'element #2 de la colonne 0: ", data.getValueDbl(2, 0))
print("Contenu de l'element #2 de la colonne 1: ", data.getValueInt(2, 1))
print("Contenu de l'element #2 de la colonne 2: ", data.getValueStr(2, 2))
data.setValueDbl(2, 0, 4)
data.setValueInt(2, 1, 8)
data.setValueStr(2, 2, "foobar")
print("Contenu de l'element #2 de la colonne 0: ", data.getValueDbl(2, 0))
print("Contenu de l'element #2 de la colonne 1: ", data.getValueInt(2, 1))
print("Contenu de l'element #2 de la colonne 2: ", data.getValueStr(2, 2))
print("Contenu du Vecteur Double:", data.getArrayDbl(0))
print("Contenu du Vecteur Int:", data.getArrayInt(1))
print("Contenu du Vecteur String:", data.getArrayStr(2))
print("Nom de la colonne 1 :", data.getName(1))
