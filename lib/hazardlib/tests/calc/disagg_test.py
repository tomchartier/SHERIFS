# The Hazard Library
# Copyright (C) 2012-2022 GEM Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import unittest
import os.path
import numpy

from openquake.hazardlib.nrml import to_python
from openquake.hazardlib.calc import disagg
from openquake.hazardlib import nrml
from openquake.hazardlib.sourceconverter import SourceConverter
from openquake.hazardlib.gsim.campbell_2003 import Campbell2003
from openquake.hazardlib.geo import Point
from openquake.hazardlib.imt import PGA, SA
from openquake.hazardlib.site import Site
from openquake.hazardlib.gsim.bradley_2013 import Bradley2013
from openquake.hazardlib import sourceconverter

DATA_PATH = os.path.dirname(__file__)


class BuildDisaggDataTestCase(unittest.TestCase):

    def test_magnitude_bins(self):
        """ Testing build disaggregation matrix """
        fname = os.path.join(DATA_PATH, 'data', 'ssm.xml')
        converter = sourceconverter.SourceConverter(50., 1., 10, 0.1, 10)
        groups = to_python(fname, converter)
        sources = []
        for g in groups:
            sources += g.sources
        site = Site(Point(172.63, -43.53), vs30=250, vs30measured=False,
                    z1pt0=330)
        imt = SA(3.0)
        iml = 0.25612220
        gsim_by_trt = {"Active Shallow Crust": Bradley2013()}
        truncation_level = 3.0
        n_epsilons = 1
        mag_bin_width = 0.1
        dist_bin_width = 100.
        coord_bin_width = 100.
        # Compute the disaggregation matrix
        edges, mtx = disagg.disaggregation(sources, site, imt, iml,
                                           gsim_by_trt, truncation_level,
                                           n_epsilons, mag_bin_width,
                                           dist_bin_width, coord_bin_width)
        tm = disagg.mag_pmf(mtx[:, :, :, :, :, 0])
        numpy.testing.assert_array_less(numpy.zeros_like(tm[2:]), tm[2:])


class DigitizeLonsTestCase(unittest.TestCase):

    def setUp(self):
        # First test
        self.lons1 = numpy.array([179.2, 179.6, 179.8, -179.9, -179.7, -179.1])
        self.bins1 = numpy.array([179.0, 179.5, 180.0, -179.5, -179])
        # Second test
        self.lons2 = numpy.array([90.0, 90.3, 90.5, 90.7, 91.3])
        self.bins2 = numpy.array([90.0, 90.5, 91.0, 91.5])

    def test1(self):
        idx = disagg._digitize_lons(self.lons1, self.bins1)
        expected = numpy.array([0, 1, 1, 2, 2, 3], dtype=int)
        numpy.testing.assert_equal(idx, expected)

    def test2(self):
        idx = disagg._digitize_lons(self.lons2, self.bins2)
        expected = numpy.array([0, 0, 1, 1, 2], dtype=int)
        numpy.testing.assert_equal(idx, expected)


class DisaggregateTestCase(unittest.TestCase):
    def setUp(self):
        d = os.path.dirname(os.path.dirname(__file__))
        source_model = os.path.join(d, 'source_model/multi-point-source.xml')
        [self.sources] = nrml.to_python(source_model, SourceConverter(
            investigation_time=50., rupture_mesh_spacing=2.))
        self.site = Site(Point(0.1, 0.1), 800, z1pt0=100., z2pt5=1.)
        self.imt = PGA()
        self.iml = 0.1
        self.truncation_level = 1
        self.trt = 'Stable Continental Crust'
        gsim = Campbell2003()
        self.gsims = {self.trt: gsim}

    def test(self):
        # a test sensitive to gsim.minimum_distance
        bin_edges, matrix = disagg.disaggregation(
            self.sources, self.site, self.imt, self.iml, self.gsims,
            self.truncation_level, n_epsilons=3,
            mag_bin_width=3, dist_bin_width=4, coord_bin_width=2.4)
        mag_bins, dist_bins, lon_bins, lat_bins, eps_bins, trt_bins = bin_edges
        aaae = numpy.testing.assert_array_almost_equal
        aaae(mag_bins, [3, 6, 9])
        aaae(dist_bins, [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52,
                         56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104,
                         108, 112])
        aaae(lon_bins, [-0.904195, 0.1, 1.104195])
        aaae(lat_bins, [-0.904194, 0.1, 1.104194])
        aaae(eps_bins, [-1, -0.3333333, 0.3333333, 1])
        self.assertEqual(trt_bins, [self.trt])
        aaae(matrix.shape, (2, 27, 2, 2, 3, 1))
        aaae(matrix.sum(), 6.14179818e-11)


class PMFExtractorsTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.aae = numpy.testing.assert_almost_equal

        # test matrix is not normalized, but that's fine for test
        self.matrix = numpy.array(
            [  # magnitude
                [  # distance
                    [  # longitude
                        [  # latitude
                            [  # epsilon
                                0.10, 0.11, 0.12],
                            [
                                0.00, 0.10, 0.20]],
                        [
                            [
                                0.74, 0.20, 0.95],
                            [
                                0.52, 0.49, 0.21]]],
                    [
                        [
                            [
                                0.16, 0.61, 0.53],
                            [
                                0.95, 0.34, 0.31]],
                        [
                            [
                                0.5, 0.61, 0.7],
                            [
                                0.40, 0.84, 0.24]]]],
                [
                    [
                        [
                            [
                                0.40, 0.32, 0.06],
                            [
                                0.47, 0.93, 0.70]],
                        [
                            [
                                0.03, 0.94, 0.12],
                            [
                                0.93, 0.13, 0.23]]],
                    [
                        [
                            [
                                0.11, 0.85, 0.85],
                            [
                                0.67, 0.84, 0.41]],
                        [
                            [
                                0.39, 0.88, 0.20],
                            [
                                0.14, 0.61, 0.67]]]]])

    def test_mag(self):
        pmf = disagg.mag_pmf(self.matrix)
        self.aae(pmf, [1.0, 1.0])

    def test_dist(self):
        pmf = disagg.dist_pmf(self.matrix)
        self.aae(pmf, [1.0, 1.0])

    def test_trt(self):
        pmf = disagg.trt_pmf(self.matrix[None])
        # NB: self.matrix.shape -> (2, 2, 2, 2, 3)
        # self.matrix[None].shape -> (1, 2, 2, 2, 2, 3)
        self.aae(pmf, [1.0])

    def test_mag_dist(self):
        pmf = disagg.mag_dist_pmf(self.matrix)
        self.aae(pmf, [[0.9989792, 0.999985], [0.9999897, 0.999996]])

    def test_mag_dist_eps(self):
        pmf = disagg.mag_dist_eps_pmf(self.matrix)
        self.aae(pmf, [[[0.88768, 0.673192, 0.972192],
                        [0.9874, 0.98393824, 0.9260596]],
                       [[0.9784078, 0.99751528, 0.8089168],
                        [0.84592498, 0.9988768, 0.976636]]])

    def test_lon_Lat(self):
        pmf = disagg.lon_lat_pmf(self.matrix)
        self.aae(pmf, [[0.9991665, 0.9999943],
                       [0.9999982, 0.9999268]])

    def test_mag_lon_lat(self):
        pmf = disagg.mag_lon_lat_pmf(self.matrix)
        self.aae(pmf, [[[0.89146822, 0.9836056],
                        [0.9993916, 0.98589012]],
                       [[0.99232001, 0.99965328],
                        [0.99700079, 0.99480979]]])

    def test_mean(self):
        # for doc purposes: the mean of PMFs is not the PMF of the mean
        numpy.random.seed(42)
        matrix = numpy.random.random(self.matrix.shape)
        pmf1 = disagg.mag_pmf(self.matrix)
        pmf2 = disagg.mag_pmf(matrix)
        mean = (matrix + self.matrix) / 2
        numpy.testing.assert_allclose(
            (pmf1 + pmf2) / 2, [1, 1])
        numpy.testing.assert_allclose(
            disagg.mag_pmf(mean), [0.99999944, 0.99999999])
