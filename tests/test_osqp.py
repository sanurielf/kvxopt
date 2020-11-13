import unittest

class TestOSQP(unittest.TestCase):

    def setUp(self):
        try:
            from kvxopt import osqp, matrix, sparse
            c = matrix([-4., -5.])
            G = sparse(matrix([[2., 1., -1., 0.], [1., 2., 0., -1.]]))
            h = matrix([3., 3., 0., 0.])
            A = sparse(matrix([1.0,1.0],(1,2)))
            b = matrix(1.0)
            self._prob_data = (c,G,h,A,b)

        except ImportError:
            self.skipTest("OSQP not available")

    def assertAlmostEqualLists(self,L1,L2,places=7, msg=None):
        self.assertEqual(len(L1),len(L2))
        for u,v in zip(L1,L2): self.assertAlmostEqual(u,v,places, msg=msg)



    def test_lp(self):
        from kvxopt import solvers, osqp
        c,G,h,A,b = self._prob_data
        sol1 = solvers.lp(c,G,h)
        self.assertTrue(sol1['status']=='optimal')
        sol2 = solvers.lp(c,G,h,solver='osqp')
        self.assertTrue(sol2['status']=='optimal')

        self.assertAlmostEqualLists(list(sol1['x']), list(sol2['x']), 2)
        self.assertAlmostEqualLists(list(sol1['z']), list(sol2['z']), 2)
        self.assertAlmostEqualLists(list(sol1['y']), list(sol2['y']), 2)

        sol3 = solvers.lp(c,G,h,A,b)
        self.assertTrue(sol3['status']=='optimal')
        sol4 = solvers.lp(c,G,h,A,b,solver='osqp')
        self.assertTrue(sol4['status']=='optimal')
        self.assertAlmostEqualLists(list(sol3['x']), list(sol4['x']), 2)
        self.assertAlmostEqualLists(list(sol3['z']), list(sol4['z']), 2)
        self.assertAlmostEqualLists(list(sol3['y']), list(sol4['y']), 2)

        sol5 = osqp.qp(c,G,h)
        self.assertTrue(sol5[0]=='solved')
        sol6 = osqp.qp(c,G,h,A,b)
        self.assertTrue(sol6[0]=='solved')
        sol7 = osqp.qp(c,G,h,None,None)
        self.assertTrue(sol7[0]=='solved')


    def test_qp(self):

        from kvxopt import solvers, osqp, spmatrix, sparse, matrix
        # Example from OSQP
        q = matrix([1., 1.])
        P = sparse(matrix([[4, 1], [0, 2]]))
        G = sparse(matrix([[1., 1, 0, -1, -1,  0],
                           [1., 0, 1, -1,  0, -1]]))
        h = matrix([1, 0.7, 0.7, -1, 0, 0])

        x = [0.3, 0.7]
        y = []
        z = [0, 0, 0.2, 2.9, 0, 0]
        obj = 1.88

        sol1 = solvers.qp(P, q, G, h, solver = 'osqp')
        self.assertTrue(sol1['status']=='optimal')

        self.assertAlmostEqualLists(list(sol1['x']), x, 2)
        self.assertAlmostEqualLists(list(sol1['y']), y, 2)
        self.assertAlmostEqualLists(list(sol1['z']), z, 1)
        self.assertAlmostEqual(sol1['primal objective'], obj, 5)


if __name__ == '__main__':
    unittest.main()
