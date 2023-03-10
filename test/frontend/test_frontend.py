import unittest
import saliweb.test

# Import the emseqfinder frontend with mocks
emseqfinder = saliweb.test.import_mocked_frontend("emseqfinder", __file__,
                                           '../../frontend')


class Tests(saliweb.test.TestCase):

    def test_index(self):
        "Test index page"
        c = emseqfinder.app.test_client()
        rv = c.get('/')
        self.assertIn('Main Page', rv.data)


if __name__ == '__main__':
    unittest.main()

