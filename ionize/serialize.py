import json
import numpy as np


def _serialize(serial, nested, compact):
    if nested:
        return serial

    if compact:
        sort_keys, indent, separators = True, None, (',', ':')
    else:
        sort_keys, indent, separators = True, 4, (', ', ': ')

    return json.dumps(serial, default=encode, sort_keys=sort_keys,
                      indent=indent, separators=separators)


def encode(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    try:
        return obj.serialize(nested=True)
    except:
        return json.JSONEncoder().default(obj)
