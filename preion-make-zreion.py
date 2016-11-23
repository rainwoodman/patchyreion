"""
    Produce a patchy reionization modulation field.
    Different parts of the universe reionizes at different times.
    Denser region is ionized earlier because there is more source.
    Thus there is a correlation between density field and the reionization
    time (time the UV background penetrates the IGM).

    This script implements a correlation model in Battaglia et al. 2013
    http://adsabs.harvard.edu/abs/2013ApJ...776...81B

    The model is used in BlueTides Simulation.

    This version of code takes the non-linear particle output of FastPM.

    It is supposed to write a UV Modulation in a format known by MP-Gadget.


    Authors:

        Yu Feng <rainwoodman@gmail.com>

"""

import argparse
import bigfile
from pmesh.pm import ParticleMesh, RealField, ComplexField
import logging
from mpi4py import MPI
import numpy

ap = argparse.ArgumentParser("preion-make-zreion.py")
ap.add_argument("fastpm", help='Non-linear particle data from FastPM')
ap.add_argument("output", help='name of bigfile to store the mesh')
ap.add_argument("--dataset", default='Zreion_table', help='name of the dataset that stores the reionization redshift')
ap.add_argument("--resolution", type=float, default=1.0, help='resolution in Mpc/h')
ap.add_argument("--filtersize", type=float, default=1.0, help='resolution in Mpc/h')
ap.add_argument("--chunksize", type=int, default=1024*1024*32, help='number of particle to read at once')
logger = logging
logging.basicConfig(level=logging.INFO)

def tophat(R, k):
    rk = R * k
    mask = rk == 0
    rk[mask] = 1
    ans = 3.0/(rk*rk*rk)*(numpy.sin(rk)-(rk)*numpy.cos(rk))
    ans[mask] = 1
    return ans

def Bk(k):
    # patchy reionization model
    # FIXME: need a citation, any meta parameters?

    b0 = 1.0 / 1.686;
    k0 = 0.185;
    al = 0.564;
    ans =  b0/pow(1 + (k/k0),al);
    return ans;

def main():
    ns = ap.parse_args()
    comm = MPI.COMM_WORLD

    ff = bigfile.BigFileMPI(comm, ns.fastpm)
    with ff['.'] as bb:
        BoxSize = bb.attrs['BoxSize'][0]
        Redshift = 1 / bb.attrs['ScalingFactor'][0] - 1

    Nmesh = int(BoxSize / ns.resolution * 2)
    # round it to 8.
    Nmesh -= Nmesh % 8

    if comm.rank == 0:
        logger.info("source = %s", ns.fastpm)
        logger.info("output = %s", ns.output)
        logger.info("BoxSize = %g", BoxSize)
        logger.info("Redshift = %g", Redshift)
        logger.info("Nmesh = %g", Nmesh)

    pm = ParticleMesh([Nmesh, Nmesh, Nmesh], BoxSize, comm=comm)

    real = RealField(pm)
    real[...] = 0

    with ff['Position'] as ds:
        logger.info(ds.size)
        for i in range(0, ds.size, ns.chunksize):
            sl = slice(i, i + ns.chunksize)
            pos = ds[sl]
            layout = pm.decompose(pos)
            lpos = layout.exchange(pos)
            real.paint(lpos, hold=True)

    mean = real.cmean()

    if comm.rank == 0:
        logger.info("mean particle per cell = %s", mean)

    real[...] /= mean
    real[...] -= 1

    complex = real.r2c()

    for k, i, slab in zip(complex.slabs.x, complex.slabs.i, complex.slabs):
        k2 = sum(kd ** 2 for kd in k)
        # tophat
        f = tophat(ns.filtersize, k2 ** 0.5)
        slab[...] *= f
        # zreion
        slab[...] *= Bk(k2 ** 0.5)
        slab[...] *= (1 + Redshift)

    real = complex.c2r()
    real[...] += Redshift

    mean = real.cmean()
    if comm.rank == 0:
        logger.info("zreion.mean = %s", mean)

    buffer = numpy.empty(real.size, real.dtype)
    real.sort(out=buffer)
    if comm.rank == 0:
        logger.info("sorted for output")

    with bigfile.BigFileMPI(comm, ns.output, create=True) as ff:
        with ff.create_from_array(ns.dataset, buffer) as bb:
            bb.attrs['BoxSize'] = BoxSize
            bb.attrs['Redshift'] = Redshift
            bb.attrs['TopHatFilterSize'] = ns.filtersize
            bb.attrs['Nmesh'] = Nmesh
        #
        # hack: compatible with current MPGadget. This is not really needed
        # we'll remove the bins later, since BoxSize and Nmesh are known.
        with ff.create("XYZ_bins", dtype='f8', size=Nmesh) as bb:
            if comm.rank == 0:
                bins = numpy.linspace(0, BoxSize * 1000., Nmesh, dtype='f8')
                bb.write(0, bins)

    if comm.rank == 0:
        logger.info("done. written at %s", ns.output)

if __name__ == '__main__':
    main()
