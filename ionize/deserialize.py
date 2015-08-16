from .Solution import Solution
from .Ion import Ion
import json

def deserialize(serial):
        serial = json.loads(serial, object_hook=object_hook)


def object_hook(obj):
    if '__ion__' in obj:
        obj.pop('__ion__')
        return Ion(**obj)
    elif '__solution__' in obj:
        obj.pop('__solution__')
        return Solution(**obj)
    return obj
