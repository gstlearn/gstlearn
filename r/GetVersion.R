################################################################################
#                                                                              #
#                            gstlearn R package                                #
#                                                                              #
# Copyright (c) (2023) MINES PARIS / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://github.com/gstlearn                                         #
# License: GPL v3                                                              #
#                                                                              #
################################################################################

r_version = paste(R.Version()$major,R.Version()$minor,sep=".")
cat(r_version)
