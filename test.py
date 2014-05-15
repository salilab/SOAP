"""
   SOAP testing module

"""

class TestMHC2(unittest.TestCase):

    def setUp(self):
        pass
        self.seq = range(10)

    def test_decoys_preparation(self):
        # make sure the shuffled sequence does not lose any elements
        random.shuffle(self.seq)
        self.seq.sort()
        self.assertEqual(self.seq, range(10))

    def test_model_selection(self):
        pass



if __name__ == '__main__':
    unittest.main()

