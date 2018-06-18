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

class DataBase(object):

    class MaskManager(object):
        def __init__(self, values):
            self._values = np.array(values)
            self._cached_masks = {}

        def get_mask(self, selector):
            """
            >>> mm = ObjectBase.MaskManager([])
            >>> list(mm.get_mask(lambda x: x != 0))
            []

            >>> list(mm.get_mask(lambda x: x == 0))
            []

            >>> mm = ObjectBase.MaskManager([1, 2])
            >>> list(mm.get_mask(lambda x: x != 0))
            [True, True]

            >>> list(mm.get_mask(lambda x: x == 0))
            [False, False]

            >>> list(mm.get_mask(lambda x: x > 1))
            [False, True]

            >>> list(mm.get_mask([]))
            [False, False]

            >>> list(mm.get_mask([1]))
            [True, False]

            >>> list(mm.get_mask([1, 2]))
            [True, True]

            >>> list(mm.get_mask([3,]))
            [False, False]

            >>> mm = ObjectBase.MaskManager([u'a', u'a', u'b'])
            >>> list(mm.get_mask(['b']))
            [False, False, True]
            """
            if hasattr(selector, '__call__'):
                return selector(self._values)
            
            return self.combine_masks(selector)

        def combine_masks(self, acceptedValues):
            if not acceptedValues:
                return np.zeros_like(self._values, dtype=bool)

            # XXX: Python3 fix
            return self._sum_masks(list(map(self._get_masks_matching_value, acceptedValues)))

        def _sum_masks(self, masks):
            if len(masks) == 1:
                return masks[0]

            return sum(masks[1:], masks[0])

        def _get_masks_matching_value(self, value):
            try:
                return self._cached_masks[value]

            except KeyError:
                return self._make_and_cache_mask(value)

        def _make_and_cache_mask(self, value):
            mask = self._values == value
            self._cachedMasks[value] = mask
            return mask

    def __init__(self, converters={}):
        """
        """
        self.__objects = np.array([], dtype=object)
        self.__cachedMaskManagers = {}
        self.__converters = dict(converters)

    def __len__(self):
        return len(self.__objects)

    def put(self, objects):
        self.__objects = np.append(self.__objects,
                                   objects if isinstance(objects, Sequence) else list(objects))
        self.__cachedMaskManagers.clear()

    def get(self, filters=None):
        return list(self.__getFilteredObjects(filters))

    def __getFilteredObjects(self, filters):
        if filters:
            return self.__objects[self.__getProductOfMasks(filters)]

        return self.__objects
      
    def __getProductOfMasks(self, selectors):
        mask = True
        for attributeName, selector in selectors.items():
            mask = mask * self.__getMask(attributeName, selector)
            
        return mask

    def __getMask(self, attributeName, selector):
        return self.__getMaskManager(attributeName).getMask(selector)

    def __getMaskManager(self, attributeName):
        try:
            return self.__cachedMaskManagers[attributeName]

        except KeyError:
            maskManager = self.MaskManager(self.__getConvertedAttributeValues(attributeName))
        self.__cachedMaskManagers[attributeName] = maskManager
        return maskManager
  
    def __getConvertedAttributeValues(self, attributeName):
        attributeValues = self.getAttributes(attributeName)
        if attributeName in self.__converters:
            # XXX: Python3 fix - makes NumPy array working
            return list(map(self.__converters[attributeName], attributeValues))

        return attributeValues

    def getAttributes(self, *attributeNames):
        """
        >>> ob = ObjectBase()
        >>> ob.put([ClassA(ClassB(1, 2), 1), ClassA(ClassB(2, 3), 2),
        ...        ClassA(ClassB(4, 3), 3)])
        >>> ob.getAttributes('a')
        [ClassB(c=1, d=2), ClassB(c=2, d=3), ClassB(c=4, d=3)]
        
        >>> ob.getAttributes('b')
        [1, 2, 3]
        
        >>> ob.getAttributes('a.c')
        [1, 2, 4]
        
        >>> ob.getAttributes('a.c', 'b')
        [(1, 1), (2, 2), (4, 3)]
        
        >>> ob = ObjectBase()
        >>> ob.getAttributes('a.c')
        []

        >>> ob = ObjectBase({'a': lambda x: x.c})
        >>> ob.put([ClassA(ClassB(1, 2), 1)])
        >>> ob.getAttributes('a')
        [ClassB(c=1, d=2)]
        """
        # XXX: Python3 fix
        return list(map(attrgetter(*attributeNames), self.__objects))
