class Animal(BaseNode):
  class DifferentMouseError(ValueError):
    pass

  __slots__ = ('Name', 'Tag', 'Sex', 'Notes')

  def __init__(self, Name, Tag, Sex=None, Notes=None):
    self.__Name = Name
    self.__Tag = Tag
    self.__Sex = Sex
    self.__Notes = Notes if isinstance(Notes, frozenset) else frozenset() if Notes is None else frozenset({unicode(Notes)})

  @classmethod
  def from_row(cls, Name, Tag, Sex=None, Notes=None):
    return cls(unicode(Name),
               frozenset({unicode(Tag.strip())}),
               None if Sex is None else unicode(Sex),
               Notes)

  def clone(self):
    return self.__class__(self.__Name,
                          self.__Tag,
                          self.__Sex,
                          self.__Notes)

  def __eq__(self, other):
    if isString(other):
      return other.__class__(self) == other

    if self.__Name != other.__Name:
      return False

    if self.__Sex == other.__Sex or other.__Sex is None or self.__Sex is None:
      return True

    return NotImplemented

  def __ne__(self, other):
    if isString(other):
      return other.__class__(self) != other

    if self.__Name != other.__Name:
      return True

    if self.__Sex == other.__Sex or other.__Sex is None or self.__Sex is None:
      return False

    return NotImplemented

  def __hash__(self):
    return self.__Name.__hash__()


  if sys.version_info >= (3, 0):
    def __str__(self):
      return self.__Name

  else:
    def __str__(self):
      return self.__Name.encode('utf-8')

    def __unicode__(self):
      return self.__Name

  def __repr__(self):
    result = '< Animal %s (' % self.Name
    if self.__Sex is not None:
      result += '%s; ' % self.Sex

    result += 'Tag: ' if len(self.__Tag) == 1 else 'Tags: '
    result += ', '.join(str(t) for t in sorted(self.__Tag))
    return result + ') >'

  def merge(self, other):
    if self != other:
      raise self.DifferentMouseError

    if self.__Sex is None:
      self.__Sex = other.__Sex

    self.__Notes = self.__Notes | other.__Notes
    self.__Tag = self.__Tag | other.__Tag
