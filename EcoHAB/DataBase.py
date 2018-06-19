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
            self.values = np.array(values)
            self.cached_masks = {}

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
                return selector(self.values)
            
            return self.combine_masks(selector)

        def combine_masks(self, acceptedValues):
            if not acceptedValues:
                return np.zeros_like(self.values, dtype=bool)

            # XXX: Python3 fix
            return self.sum_masks(list(map(self.get_cached_masks, acceptedValues)))

        def sum_masks(self, masks):
            if len(masks) == 1:
                return masks[0]

            return sum(masks[1:], masks[0])

        def get_cached_masks(self, value):
            try:
                return self.cached_masks[value]

            except KeyError:
                return self.make_and_cache_mask(value)

        def make_and_cache_mask(self, value):
            mask = self.values == value
            self.cachedMasks[value] = mask
            return mask

    def __init__(self, converters={}):
        """
        """
        self.objects = np.array([], dtype=object)
        self.cachedMaskManagers = {}
        self.converters = dict(converters)

    def __len__(self):
        return len(self.objects)

    def put(self, objects):
        self.objects = np.append(self.objects,
                                   objects if isinstance(objects, Sequence) else list(objects))
        self.cachedMaskManagers.clear()

    def get(self, filters=None):
        return list(self.get_filtered_objects(filters))

    def get_filtered_objects(self, filters):
        if filters:
            return self.objects[self.apply_mask(filters)]

        return self.objects
      
    def apply_mask(self, selectors):
        mask = True
        for attributeName, selector in selectors.items():
            mask = mask * self.get_mask(attributeName, selector)
            
        return mask

    def get_mask(self, attributeName, selector):
        return self.get_mask_manager(attributeName).get_mask(selector)

    def get_mask_manager(self, attributeName):
        try:
            return self.cached_mask_managers[attributeName]

        except KeyError:
            maskManager = self.MaskManager(self.get_converted_attribute_values(attributeName))
        self.cached_mask_managers[attributeName] = maskManager
        return maskManager
  
    def get_converted_attribute_values(self, attributeName):
        attributeValues = self.get_attributes(attributeName)
        if attributeName in self.converters:
            # XXX: Python3 fix - makes NumPy array working
            return list(map(self.converters[attributeName], attributeValues))

        return attributeValues

    def get_attributes(self, *attributeNames):
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
        return list(map(attrgetter(*attributeNames), self.objects))
