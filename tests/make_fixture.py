"""Reproduce the included ITK HDF5 fixture. Requires numpy and h5py only here,
not for building or running the R/C++ implementation.
The composite applies a constant +1 mm LPS-x displacement, then +1 mm affine.
"""
from pathlib import Path
import h5py
import numpy as np
path = Path(__file__).parent / 'fixtures/composite_translation.h5'
path.parent.mkdir(exist_ok=True)
with h5py.File(path, 'w') as f:
    string = h5py.string_dtype('ascii')
    for key, value in [('HDFVersion', h5py.version.hdf5_version),
                       ('ITKVersion', '5.4.0'), ('OSName', 'fixture'), ('OSVersion', '1')]:
        f.create_dataset(key, data=[value], dtype=string)
    group = f.create_group('TransformGroup')
    g = group.create_group('0')
    g.create_dataset('TransformType', data=['CompositeTransform_double_3_3'], dtype=string)
    g = group.create_group('1')
    g.create_dataset('TransformType', data=['AffineTransform_double_3_3'], dtype=string)
    g.create_dataset('TransformFixedParameters', data=np.zeros(3))
    g.create_dataset('TransformParameters', data=np.r_[np.eye(3).ravel(), [1.,0.,0.]])
    g = group.create_group('2')
    g.create_dataset('TransformType', data=['DisplacementFieldTransform_double_3_3'], dtype=string)
    # Field covers the entire input/output domain, including repeated transformations.
    dims = [32,32,32]
    fixed = np.r_[dims, [-16.,-16.,-16.], [1.,1.,1.], np.eye(3).ravel()]
    g.create_dataset('TransformFixedParameters', data=fixed)
    displacement = np.tile([1.,0.,0.], np.prod(dims))
    g.create_dataset('TransformParameters', data=displacement)
print(path)
