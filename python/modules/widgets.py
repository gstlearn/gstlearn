################################################################################
#                                                                              #
#                         gstlearn Python package                              #
#                                                                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://gstlearn.org                                                #
# License: BSD 3-clause                                                        #
#                                                                              #
################################################################################
import gstlearn as gl
import matplotlib.pyplot as plt
import ipywidgets as ipw
import numpy as np
import re
from IPython.display import display
from traitlets import Unicode

def boolean(title="", value=False, eventhandler=None):
    object = ipw.Checkbox(value=value, description=title, disabled=False, indent=False)
    if eventhandler is not None:
        object.observe(eventhandler, names='value')

    return object

def sliderInt(title="", value=0, mini=0, maxi=100, step=1, eventhandler=None):
    object  = ipw.IntSlider(min=mini, max=maxi, step=step, value=value,
                                description=title, continuous_update=False)
    if eventhandler is not None:
        object.observe(eventhandler, names='value')
    return object
        
def sliderFloat(title="", value=0., mini=0., maxi=100., step=1., eventhandler=None):
    object  = ipw.FloatSlider(min=mini, max=maxi, step=step, value=value,
                                description=title, continuous_update=False)
    if eventhandler is not None:
        object.observe(eventhandler, names='value')
    return object

def text(title="", value='', placeholder=''):
    object = ipw.Text(value=value, placeholder=placeholder,description=title,
                      disabled=False)
    return object

def label(title="", value=''):
    object = ipw.Label(value=value, description=title)
    return object

def dropDown(title="", options=['1', '2', '3'], value=0, eventhandler=None):
    object = ipw.Dropdown(options = options, value = value, description=title,  
                          continuous_update=False)
    if eventhandler is not None:
        object.observe(eventhandler, names='value')
    return object

class WModel(ipw.VBox):
    value = Unicode()

    def __init__(self, models, changeCallback, defrank=0, defrange=20, defsill=1, defparam=1, **kwargs):
        defvalue = list(models.keys())[defrank]
        self.dropType      = dropDown(title='Type', options = models.keys(), value=defvalue,
                                      eventhandler=self.__update_value)
        self.sliderRangeX  = sliderInt(title='RangeX', value=defrange, mini=1, maxi=50,
                                       eventhandler=self.__update_value)
        self.sliderRangeY  = sliderInt(title='RangeY', value=defrange, mini=1, maxi=50,
                                       eventhandler=self.__update_value)
        self.sliderParam   = sliderInt(title='Parameter', value=defparam, mini=0, maxi=4,
                                       eventhandler=self.__update_value)
        self.sliderAngle   = sliderInt(title='Angle', value=0, mini=0, maxi=180,
                                       eventhandler=self.__update_value)
        self.sliderSill    = sliderInt(title='Sill', value=defsill, mini=0, maxi=100,
                                       eventhandler=self.__update_value)
        
        self.hb1 = ipw.HBox([self.sliderRangeX, 
                             self.sliderRangeY, 
                             self.sliderParam])
        self.hb2 = ipw.HBox((self.sliderAngle, self.sliderSill))
        self.changeCallback = changeCallback
        self.models = models
        self.monmodel = gl.Model.createFromParam(type=models[defvalue],range=defrange,sill=defsill,
                                                 param=defparam)

#        self.__update_value()
        self.observe(self.__update_children, names='value')

        super().__init__(children=[self.dropType, 
                                   self.hb1, self.hb2],
                         **kwargs)

    def __toModel(self):
        self.monmodel = gl.Model.createFromParam(type=self.models[self.dropType.value],
                                                 ranges=[self.sliderRangeX.value,self.sliderRangeY.value],
                                                 angles=[self.sliderAngle.value,0],
                                                 param=self.sliderParam.value,
                                                 sill=self.sliderSill.value)
        
    def __toValue(self, myModel):
        cova   = myModel.getCovAniso(0)
        type   = cova.getType()
        rangeX = cova.getRange(0)
        rangeY = cova.getRange(1)
        param  = cova.getParam()
        angle  = cova.getAnisoAngles()[0]
        sill   = cova.getSill(0,0)
        self.value = "{},{},{},{},{},{}".format(type, rangeX, rangeY, param, angle, sill)
        
    def __update_children(self, *args):
        numvals = self.value.split(',')
        self.dropType.value         = numvals[0]
        self.sliderRangeX.value     = float(numvals[1])
        self.sliderRangeY.value     = float(numvals[2])
        self.sliderParam.value      = float(numvals[3])
        self.sliderAngle.value      = float(numvals[4])
        self.sliderSill.value       = float(numvals[5])
        self.__toModel()

    def __update_value(self, *args):
        self.value = "{},{},{},{},{},{}".format(self.dropType.value,
                                                self.sliderRangeX.value, 
                                                self.sliderRangeY.value,
                                                self.sliderParam.value,
                                                self.sliderAngle.value,
                                                self.sliderSill.value)
        self.__toModel()
        self.changeCallback(self.monmodel)

    def update_model(self, myModel):
        self.__toValue(myModel)
