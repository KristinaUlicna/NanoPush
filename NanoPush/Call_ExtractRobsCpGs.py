from ExtractUnMethCpGs_Function import ExtractCpGsNewModel

old_unmeth_U, old_meth_U = ExtractCpGsNewModel("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/OldModel/unmeth/")
print ("\nUnmeth_U:", len(old_unmeth_U), "\n\tMeth_U:", len(old_meth_U))

old_unmeth_M, old_meth_M = ExtractCpGsNewModel("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/OldModel/meth/")
print ("\nUnmeth_M:", len(old_unmeth_M), "\n\tMeth_M:", len(old_meth_M))


new_unmeth_U, new_meth_U = ExtractCpGsNewModel("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewModel/unmeth/")
print ("\nUnmeth_U:", len(new_unmeth_U), "\n\tMeth_U:", len(new_meth_U))

new_unmeth_M, new_meth_M = ExtractCpGsNewModel("/Users/kristinaulicna/Documents/Rotation_1/Archive/DRONE/NewModel/meth/")
print ("\nUnmeth_M:", len(new_unmeth_M), "\n\tMeth_M:", len(new_meth_M))


