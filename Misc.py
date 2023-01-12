
import numpy as np
from scipy.interpolate import interp1d, UnivariateSpline  

    
class LinearExtra(interp1d): # linear interpolation with flat extension
    def __init__(self, x, y):
        return super().__init__(x, y, kind='linear', bounds_error=False, fill_value='extrapolate')

class LinearFlat(interp1d): # linear interpolation with flat extension
    def __init__(self, x, y):
        return super().__init__(x, y, kind='linear', bounds_error=False, fill_value=(y[0], y[-1]))


class HashableArray(np.ndarray):
    def __new__(cls, *args):
        return np.asarray(*args).view(cls)     

    def __eq__(self, other):
        return np.array_equal(self, other)

    def __hash__(self):
        return hash(self.tostring())

    def __getitem__(self, *args):
        return HashableArray(super().__getitem__(*args))    


    @classmethod
    def fetch(cls, id): # id as string, e.g. 'US' / 'TGT'
        return cls.__dict__[id]   


def time_this(fn):
    """
    decorator that reports the execution time.
    """
    # use @wraps(fn) to expose the function name and docstring to the caller of
    # the decorated function (otherwise, we would see the function name and
    # docstring for the decorator, not the function it decorates).
    @wraps(fn) 
    def wrapper(*args, **kwargs):
        print('Start ' + fn.__name__)
        start = time.time()
        result = fn(*args, **kwargs)
        print(fn.__name__, 'finished in', '{0:.2f}'.format(time.time() - start), 'seconds\n')
        return result
    return wrapper 



if __name__ == '__main__':
    pass
