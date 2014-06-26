"""
   SOAP error/exception types.

"""

class NetworkError(Exception):
    """
    Network error
    """
    pass

class NanInScore(Exception):
    """
    Nan error
    """
    pass

class Bugs(Exception):
    """
    Bug in code
    """
    pass

class FatalError(Exception):
    """
    Fatal error
    """
    pass

class OutOfRange(Exception):
    """
    Index out of range error
    """
    pass
