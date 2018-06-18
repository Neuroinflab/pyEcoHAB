#!/usr/bin/env python
# encoding: utf-8
###############################################################################
#                                                                             #
#    EcoHAB library                                                           #
#                                                                             #
#    Copyright (C) 2018 Joanna Jędrzejewska-Szmek, S. Łęski, Jan Mąka         #
#    Jakub M. Dzik a.k.a. Kowalski,                                           #
#    (Laboratory of Neuroinformatics; Nencki Institute of Experimental        #
#    Biology of Polish Academy of Sciences)                                   #
#                                                                             #
#    This software is free software: you can redistribute it and/or modify    #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This software is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this software.  If not, see http://www.gnu.org/licenses/.     #
#                                                                             #
###############################################################################
from __future__ import division, print_function, absolute_import
import numpy as np
from operator import attrgetter
from collections import Sequence

class __DataBase(object):

    class MaskManager(object):
        def __init__(self, values):
            self.__values = np.array(values)
            self.__cachedMasks = {}

        def getMask(self, selector):
            """
            >>> mm = ObjectBase.MaskManager([])
            >>> list(mm.getMask(lambda x: x != 0))
            []

            >>> list(mm.getMask(lambda x: x == 0))
            []

            >>> mm = ObjectBase.MaskManager([1, 2])
            >>> list(mm.getMask(lambda x: x != 0))
            [True, True]

            >>> list(mm.getMask(lambda x: x == 0))
            [False, False]

            >>> list(mm.getMask(lambda x: x > 1))
            [False, True]

            >>> list(mm.getMask([]))
            [False, False]

            >>> list(mm.getMask([1]))
            [True, False]

            >>> list(mm.getMask([1, 2]))
            [True, True]

            >>> list(mm.getMask([3,]))
            [False, False]

            >>> mm = ObjectBase.MaskManager([u'a', u'a', u'b'])
            >>> list(mm.getMask(['b']))
            [False, False, True]
            """
            if hasattr(selector, '__call__'):
                return selector(self.__values)
            
            return self.__combineMasks(selector)

        def __combineMasks(self, acceptedValues):
            if not acceptedValues:
                return np.zeros_like(self.__values, dtype=bool)

            # XXX: Python3 fix
            return self.__sumMasks(list(map(self.__getMasksMatchingValue, acceptedValues)))

        def __sumMasks(self, masks):
            if len(masks) == 1:
                return masks[0]

            return sum(masks[1:], masks[0])

        def __getMasksMatchingValue(self, value):
            try:
                return self.__cachedMasks[value]

            except KeyError:
                return self.__makeAndCacheMask(value)

        def __makeAndCacheMask(self, value):
            mask = self.__values == value
            self.__cachedMasks[value] = mask
            return mask

