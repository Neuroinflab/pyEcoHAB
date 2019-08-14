from __future__ import division, print_function
from operator import attrgetter, methodcaller
import sys

def makePrivateSlots(attributes, name):
  prefix = '_%s__' % name
  return tuple(prefix + s for s in attributes)


class BaseNodeMetaclass(type):
  def __new__(mcl, name, bases, attrs):
    attributes = attrs['__slots__']
    slots = makePrivateSlots(attributes, name)
    attrs['__slots__'] = slots
    attrs.update(zip(attributes,
                     (property(attrgetter(s)) for s in slots)))

    return type.__new__(mcl, name, bases, attrs)
  

def BaseNode_del_(self):
  for cls in self.__class__.__mro__:
    if not hasattr(cls, '__slots__'):
      continue

    for attr in cls.__slots__:
      try:
        delattr(self, attr)

      except AttributeError:
        pass

class BaseNode(object, metaclass=BaseNodeMetaclass):
  __slots__ = ()
  _del_ = BaseNode_del_
  
    
class AntennaReadOut(object):
    __slots__ = ('EventId', 'AntennaId', 'Animal', 'Start', 'Date', 'Duration', 'source')
    def __init__(self, EventId, AntennaId, Animal, Start, Date, Duration):
        self.EventId = EventId
        self.AntennaId = AntennaId
        self.Animal = Animal
        self.Start = Start
        self.Date = Date
        self.Duration = Duration
 
    def clone(self, animalManager):
        #animal = animalManager[self.Animal]
        return self.__class__(self.EventId,
                              self.AntennaId,
                              self.Animal,
                              self.Start,
                              self.Date,
                              self.Duration)
    
    def _del_(self):
        super(Visit, self)._del_()

    def __repr__(self):
        return '<Antenna Readout %d of animal %s at %s>' % (self.AntennaId,
                                                            self.Animal,
                                                            self.Date)
class Animal(object):
  class DifferentMouseError(ValueError):
    pass

  __slots__ = ('Tag')

  def __init__(self, Tag):
    self.Tag = Tag

  @classmethod
  def fromRow(cls, Tag):
    return cls(frozenset({unicode(Tag.strip())}))

  def clone(self):
      return self.__class__(self.Tag)

  def __eq__(self, other):
    if isString(other):
      return other.__class__(self) == other

    if self.__Tag != other._Tag:
      return False

  def __ne__(self, other):
    if isString(other):
      return other.__class__(self) != other

    if self.Tag != other.Tag:
      return True

    return NotImplemented

  def __hash__(self):
    return self.Tag.__hash__()


  if sys.version_info >= (3, 0):
    def __str__(self):
      return self.Tag

  else:
    def __str__(self):
      return self.Tag.encode('utf-8')

    def __unicode__(self):
      return self.Tag

  def __repr__(self):
    result = '< Animal %s >' % self.Tag
    return result

  def merge(self, other):
    if self != other:
      raise self.DifferentMouseError
    self.Tag = self.Tag | other.Tag

